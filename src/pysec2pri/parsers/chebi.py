"""ChEBI SDF file parser for secondary-to-primary identifier mappings."""

from __future__ import annotations

from pathlib import Path

from tqdm import tqdm

from pysec2pri.models import (
    ChEBIMapping,
    MappingCardinality,
    MappingSet,
    compute_cardinality,
)
from pysec2pri.parsers.base import BaseParser


class ChEBIParser(BaseParser):
    """Parser for ChEBI SDF files.

    Extracts secondary-to-primary ChEBI identifier mappings and
    name-to-synonym relationships from ChEBI SDF format files.

    Note: SDF is not a tabular format, so Polars is not applicable here.
    """

    datasource_name = "ChEBI"
    default_source_url = (
        "https://ftp.ebi.ac.uk/pub/databases/chebi/SDF/ChEBI_complete_3star.sdf"
    )

    def parse(self, input_path: Path | str) -> MappingSet:
        """Parse a ChEBI SDF file.

        Args:
            input_path: Path to the ChEBI SDF file
                        (e.g., ChEBI_complete_3star.sdf).

        Returns:
            MappingSet with ChEBI identifier and synonym mappings.
        """
        input_path = Path(input_path)

        # Collect raw mappings first
        raw_id_mappings: list[tuple[str, str]] = []  # (primary, secondary)
        raw_name_mappings: list[tuple[str, str, str]] = []  # (primary, name, syn)

        # Read and parse the SDF file
        with input_path.open("r", encoding="utf-8") as f:
            lines = f.readlines()

        current_primary_id: str | None = None
        current_name: str | None = None
        i = 0
        total_lines = len(lines)

        # Use progress bar for line processing
        pbar = None
        if self.show_progress:
            pbar = tqdm(total=total_lines, desc="Parsing ChEBI SDF")

        while i < total_lines:
            line = lines[i].strip()

            if pbar:
                pbar.update(1)

            # Extract primary ID
            if line.startswith("> <ChEBI ID>"):
                i += 1
                if pbar:
                    pbar.update(1)
                if i < total_lines:
                    current_primary_id = lines[i].strip()

            # Extract metabolite name
            elif line.startswith("> <ChEBI Name>"):
                i += 1
                if pbar:
                    pbar.update(1)
                if i < total_lines:
                    current_name = lines[i].strip()

            # Extract secondary IDs
            elif line.startswith("> <Secondary ChEBI ID>"):
                i += 1
                if pbar:
                    pbar.update(1)
                while i < total_lines:
                    sec_line = lines[i].strip()
                    if sec_line.startswith("CHEBI:"):
                        if current_primary_id:
                            raw_id_mappings.append((current_primary_id, sec_line))
                        i += 1
                        if pbar:
                            pbar.update(1)
                    else:
                        break
                continue  # Don't increment i again

            # Extract synonyms
            elif line.startswith("> <Synonyms>"):
                i += 1
                if pbar:
                    pbar.update(1)
                while i < total_lines:
                    syn_line = lines[i].strip()
                    if syn_line and not syn_line.startswith(">"):
                        if current_primary_id and current_name:
                            raw_name_mappings.append(
                                (current_primary_id, current_name, syn_line)
                            )
                        i += 1
                        if pbar:
                            pbar.update(1)
                    else:
                        break
                continue

            # End of record
            elif line.startswith("$$$$"):
                current_primary_id = None
                current_name = None

            i += 1

        if pbar:
            pbar.close()

        # Compute cardinality for ID mappings
        cardinality_map = compute_cardinality(raw_id_mappings)

        # Create mapping objects
        mappings: list[ChEBIMapping] = []

        # Add ID mappings
        id_iter = raw_id_mappings
        if self.show_progress:
            id_iter = tqdm(raw_id_mappings, desc="Creating ID mappings")

        for primary_id, secondary_id in id_iter:
            cardinality = cardinality_map.get(
                (primary_id, secondary_id),
                MappingCardinality.ONE_TO_ONE,
            )
            mapping = ChEBIMapping(
                primary_id=primary_id,
                secondary_id=secondary_id,
                mapping_cardinality=cardinality,
                source_url=self.default_source_url,
            )
            mappings.append(mapping)

        # Add name/synonym mappings (no cardinality for synonyms)
        name_iter = raw_name_mappings
        if self.show_progress:
            name_iter = tqdm(raw_name_mappings, desc="Creating synonym mappings")

        for primary_id, name, synonym in name_iter:
            mapping = ChEBIMapping(
                primary_id=primary_id,
                secondary_id=None,
                primary_label=name,
                secondary_label=synonym,
                source_url=self.default_source_url,
            )
            mappings.append(mapping)

        # Build the mapping set
        mapping_set = MappingSet(
            mappings=mappings,
            datasource_name=self.datasource_name,
            version=self.version,
            mapping_set_id="omicsfixid_chebi_01",
            mapping_set_description=(
                "Secondary to primary ID mappings for ChEBI database, "
                "generated for the omicsFixID project."
            ),
            comment=(
                "object_label represents a synonym name by which "
                "the subject_label is also known."
            ),
            curie_map={
                "CHEBI:": "http://purl.obolibrary.org/obo/CHEBI_",
                "IAO:": "http://purl.obolibrary.org/obo/IAO_",
                "oboInOwl:": "http://www.geneontology.org/formats/oboInOwl#",
            },
        )

        return mapping_set


__all__ = ["ChEBIParser"]
