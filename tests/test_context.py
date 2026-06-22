"""Edge-case tests for pysec2pri.context (label/id/xref disambiguation)."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest
from sssom_schema import Mapping

from pysec2pri.api import crosswalk
from pysec2pri.context import ContextSpec, XrefMapping, XrefRecord, resolve_ambiguous_with_xref
from pysec2pri.parsers.base import LabelMappingSet
from pysec2pri.update_ids import update_labels

_MJ = "semapv:BackgroundKnowledgeBasedMatching"


@pytest.fixture
def label_ms() -> LabelMappingSet:
    """Build a tiny label mapping set with one ambiguous symbol, SPN/DEAF1-shaped.

    ``X`` is simultaneously the *current* symbol of gene HGNC:1 and a
    *previous/alias* symbol of gene HGNC:3 (current symbol ``Y``).
    """
    mappings = [
        Mapping(
            subject_id="HGNC:1",
            subject_label="X",
            object_id="HGNC:1",
            object_label="X",
            predicate_id="skos:exactMatch",
            mapping_justification=_MJ,
        ),
        Mapping(
            subject_id="HGNC:2",
            subject_label="X",
            object_id="HGNC:3",
            object_label="Y",
            predicate_id="IAO:0100001",
            mapping_justification=_MJ,
        ),
        Mapping(
            subject_id="HGNC:4",
            subject_label="OLD_Z",
            object_id="HGNC:5",
            object_label="Z",
            predicate_id="IAO:0100001",
            mapping_justification=_MJ,
        ),
    ]
    ms = LabelMappingSet(
        mappings=mappings,
        mapping_set_id="test:ambiguous",
        license="https://creativecommons.org/publicdomain/zero/1.0/",
    )
    object.__setattr__(ms, "_primary_labels", {"X": {"HGNC:1"}, "Y": {"HGNC:3"}, "Z": {"HGNC:5"}})
    return ms


@pytest.mark.parametrize(
    ("xref_object_id", "xref_object_label", "expected"),
    [
        ("HGNC:3", "Y", ("Y", "HGNC:3")),  # xref points to the secondary's target
        ("HGNC:1", "X", ("X", "HGNC:1")),  # xref points to the token's own identity
        ("HGNC:999", "Q", ("", None)),  # xref points to neither: a third entity
    ],
)
def test_xref_resolution_decision_matrix(
    xref_object_id: str, xref_object_label: str, expected: tuple[str, str | None]
) -> None:
    """The three possible xref outcomes for an ambiguous token: target, own identity, or neither."""
    lkp = {"X": "Y"}  # secondary 'X' -> primary 'Y'
    token_to_id = {"X": "HGNC:1", "Y": "HGNC:3"}
    record = XrefRecord(
        subject_id="ENSG1", object_id=xref_object_id, object_label=xref_object_label
    )
    index = {"ENSG1": [record]}
    token, tid, decision = resolve_ambiguous_with_xref("X", "ENSG1", lkp, index, token_to_id)
    assert (token, tid) == expected
    assert decision.accepted == (expected != ("", None))


def test_confident_cells_never_touched_despite_misleading_xref(label_ms: LabelMappingSet) -> None:
    """A non-ambiguous symbol resolves without consulting xref evidence at all."""
    df = pd.DataFrame({"ensembl": ["ENSG_ANYTHING"], "gene_name": ["OLD_Z"]})
    # If this xref were consulted, it would point the row at the *wrong* gene.
    xref_mapping = XrefMapping(
        records=[XrefRecord(subject_id="ENSG_ANYTHING", object_id="HGNC:1", object_label="X")]
    )
    out = update_labels(df, label_ms, at="gene_name", xref="ensembl", xref_mapping=xref_mapping)
    assert out["gene_name_current"].iloc[0] == "Z"


def test_context_specs_fall_through_on_failure(label_ms: LabelMappingSet, tmp_path: Path) -> None:
    """An unhelpful symbol hint falls through to a successful xref attempt; both are logged."""
    df = pd.DataFrame(
        {"gene_name": ["X"], "hint": ["NOT_A_REAL_ALIAS"], "ensembl": ["ENSG_TARGET"]}
    )
    xref_mapping = XrefMapping(
        records=[XrefRecord(subject_id="ENSG_TARGET", object_id="HGNC:3", object_label="Y")]
    )
    specs = [
        ContextSpec(kind="label", column="hint"),
        ContextSpec(kind="xref", column="ensembl", xref_mapping=xref_mapping),
    ]
    report_path = tmp_path / "multi_context.tsv"
    out = update_labels(df, label_ms, at="gene_name", context=specs, report_path=report_path)
    assert out["gene_name_current"].iloc[0] == "Y"
    log = pd.read_csv(report_path, sep="\t")
    assert list(log["stage"]) == ["label", "xref_filter"]
    assert list(log["accepted"]) == [False, True]


def test_crosswalk_ambiguous_target_left_blank(tmp_path: Path) -> None:
    """A crosswalk token with two distinct targets is ambiguous, not guessed."""
    xref_mapping = XrefMapping(
        records=[
            XrefRecord(subject_id="ENSG_DUP", object_id="HGNC:1", object_label="X"),
            XrefRecord(subject_id="ENSG_DUP", object_id="HGNC:2", object_label="W"),
        ]
    )
    report_path = tmp_path / "dup.tsv"
    result = crosswalk(
        "ENSG_DUP", frm="ensembl", to="hgnc_id", xref_mapping=xref_mapping, report_path=report_path
    )
    assert result == {"ENSG_DUP": ""}
    log = pd.read_csv(report_path, sep="\t")
    assert "ambiguous crosswalk" in log["reason"].iloc[0]
