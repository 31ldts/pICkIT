"""Backward-compatibility shim.

Historically ``AnalyzeInteractions``, ``InteractionData`` and every custom
exception lived directly in this file. After the Fase 2 modularization
(see plan-accion-modularizacion-pickit.md) they live in dedicated modules
(``analyzer.py``, ``models.py``, ``exceptions.py``, ``constants.py``, and
the ``*_mixin.py`` files). This module re-exports all of them under their
original names and original import path so that any existing code doing
``from pickit.analyze_interactions import AnalyzeInteractions`` (as
tests/conftest.py and every Fase 0 baseline test does) keeps working
unchanged — this is the plan's own exit criterion for Fase 2: "AnalyzeInteractions
sigue siendo importable igual que antes (from pickit.analyze_interactions
import AnalyzeInteractions o el nuevo path equivalente)".
"""

from .analyzer import AnalyzeInteractions
from .constants import (
    AMINO_ACID_CODES,
    ARPEGGIO_COLORS,
    ARPEGGIO_CONT,
    ARPEGGIO_INT_ENT,
    ARPEGGIO_TYPE,
    COLORS,
    DEFAULT_TEMPLATE_FILE,
    DIFF_DELIM,
    EMPTY_CELL,
    EMPTY_DASH_CELL,
    GROUP_DELIM,
    HEATMAP_MODES,
    INTERACTION_LABELS,
    PROGRAM_MODES,
    SAME_DELIM,
    is_not_empty_or_dash,
)
from .exceptions import (
    FileOrDirectoryException,
    HeatmapActivityException,
    InvalidAxisException,
    InvalidColorException,
    InvalidFileExtensionException,
    InvalidFilenameException,
    InvalidModeException,
    MissedActivityException,
    TypeMismatchException,
)
from .models import InteractionData

__all__ = [
    "AnalyzeInteractions",
    "InteractionData",
    "INTERACTION_LABELS",
    "COLORS",
    "PROGRAM_MODES",
    "HEATMAP_MODES",
    "ARPEGGIO_INT_ENT",
    "ARPEGGIO_CONT",
    "ARPEGGIO_TYPE",
    "ARPEGGIO_COLORS",
    "SAME_DELIM",
    "DIFF_DELIM",
    "GROUP_DELIM",
    "EMPTY_CELL",
    "EMPTY_DASH_CELL",
    "AMINO_ACID_CODES",
    "DEFAULT_TEMPLATE_FILE",
    "is_not_empty_or_dash",
    "TypeMismatchException",
    "FileOrDirectoryException",
    "InvalidColorException",
    "InvalidFileExtensionException",
    "InvalidAxisException",
    "InvalidFilenameException",
    "HeatmapActivityException",
    "InvalidModeException",
    "MissedActivityException",
]
