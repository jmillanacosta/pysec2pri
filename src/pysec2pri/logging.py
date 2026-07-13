"""Logging configuration for pysec2pri.

Provides the ``pysec2pri`` package logger, built by :mod:`mapkgsutils.logging`.
Only CRITICAL messages show by default; call :func:`set_log_level` to raise
verbosity.
"""

from __future__ import annotations

from mapkgsutils.logging import LOG_LEVELS, get_logger
from mapkgsutils.logging import set_log_level as _set_log_level

__all__ = [
    "LOG_LEVELS",
    "logger",
    "set_log_level",
]

logger = get_logger("pysec2pri")


def set_log_level(level: str | int) -> None:
    """Set the pysec2pri log level.

    Args:
        level: Level name (``"debug"``/``"info"``/``"warning"``/``"error"``/
            ``"critical"``) or an integer such as ``logging.INFO``.
    """
    _set_log_level(level, logger)
