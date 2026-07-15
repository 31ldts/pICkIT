"""Final assembly of ``AnalyzeInteractions`` as a composition of mixins.

Fase 2, paso 8 del plan de modularización — el último paso de la
modularización. Todo el código funcional ya vive en su módulo dedicado
(``_validation.py``, ``filter_mixin.py``, ``export_mixin.py``,
``plot_mixin.py``, ``io_mixin.py``); este fichero solo define la clase que
los compone y los atributos de configuración compartidos (constantes de
clase, ``__init__``).

Verified per the plan's exit criterion: ``dir(AnalyzeInteractions)`` before
this modularization (55 baseline tests passing against the single-file,
1376-statement original) and after (103 tests passing against 8 modules)
expose exactly the same public methods — see
plan-accion-modularizacion-pickit.md for the before/after comparison log.
"""

import os

from ._validation import ValidationMixin
from .constants import AMINO_ACID_CODES, COLORS, INTERACTION_LABELS
from .export_mixin import ExportMixin
from .filter_mixin import FilterMixin
from .io_mixin import IOMixin
from .plot_mixin import PlotMixin


class AnalyzeInteractions(ValidationMixin, FilterMixin, ExportMixin, PlotMixin, IOMixin):
    """Analyze, filter, export and plot protein-ligand interaction matrices.

    Composed from five mixins, each covering one concern:

    - ``ValidationMixin``: shared input-validation helpers.
    - ``FilterMixin``: filtering, sorting and reshaping interaction matrices.
    - ``ExportMixin``: DataFrame/Excel export.
    - ``PlotMixin``: heatmaps, bar charts and pie charts.
    - ``IOMixin``: directory/config management and Arpeggio/IChem file parsing.
    """

    ROWS = "rows"
    COLUMNS = "columns"

    ARPEGGIO = "arpeggio"
    ICHEM = "ichem"

    INPUT = "input"
    OUTPUT = "output"

    MINIMUM = "min"
    MAXIMUM = "max"
    MEAN = "mean"
    COUNT = "count"
    PERCENT = "percent"

    MAIN = "<main>"
    SIDE = "<side>"

    LOWER = "lower"
    UPPER = "upper"

    INFO = "info"
    WARNING = "warning"
    ERROR = "error"

    def __init__(self):
        """Initializes the class with the current working directory and default settings."""
        self.saving_directory = os.getcwd()  # Set the default saving directory
        self.input_directory = os.getcwd()  # Set the default input directory
        self.interaction_labels = INTERACTION_LABELS  # Default interaction labels
        self.codes = True
        self.plot_colors = COLORS  # Default color configuration
        self.plot_max_cols = 80
        self.aa = AMINO_ACID_CODES
        self.heat_max_cols = 30
        self.heat_colors = "RdYlGn"
