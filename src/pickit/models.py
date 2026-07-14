"""Data model for pickit interaction matrices.

Extracted verbatim from the original monolithic ``analyze_interactions.py``
(Fase 2, paso 2 del plan de modularización). No logic changes were made
during this extraction — only the physical location of the code changed.

Decisión de estructura de datos (Anexo B del plan): se mantiene la matriz
como lista de listas de strings (``list[list[str]]``) dentro de
``InteractionData.matrix``. Ver ``docs/decisions/data-structure.md``
(pendiente de crear cuando se extraiga ``filter_mixin.py``, que es donde se
ejercen de verdad las operaciones — transponer, indexar por residuo,
filtrar — sobre las que se apoya esa decisión) para el razonamiento
completo; por ahora esto es solo el traslado de la clase, no un rediseño.
"""


class InteractionData:
    """Container for a processed interaction matrix and its metadata.

    Attributes:
        colors (list[str]): Hex colors associated with each interaction label.
        interactions (list[str]): Interaction type labels, in the same order
            used to encode interaction codes inside matrix cells.
        ligand (bool): Whether ligand atoms were included when the matrix
            was built.
        matrix (list[list[str]]): The interaction matrix itself, header row
            + labeled rows/columns, as produced by ``analyze_files`` and
            consumed by the filter/export/plot methods.
        mode (str): The program mode used to build the matrix (e.g.
            ``"arpeggio"`` or ``"ichem"``).
        protein (bool): Whether protein atoms were included when the matrix
            was built.
        subunit (bool): Whether protein subunits were kept distinct when the
            matrix was built.
        subunits_set (set[str]): The set of subunit identifiers present in
            the matrix.
    """

    def __init__(self, colors, interactions, ligand, matrix, mode, protein, subunit, subunits_set):
        self.colors = colors
        self.interactions = interactions
        self.ligand = ligand
        self.matrix = matrix
        self.mode = mode
        self.protein = protein
        self.subunit = subunit
        self.subunits_set = subunits_set

    def compare(self, other):
        """Compare this InteractionData against another object attribute by attribute.

        Args:
            other: Any object. Doesn't need to be an InteractionData instance,
                only to expose the same attributes (duck typing, matching the
                original implementation's intent).

        Returns:
            str | dict: ``"There are no changes."`` if every attribute
            matches, a dict of ``{attribute: (self_value, other_value)}``
            for attributes that differ, or an explanatory string if
            ``other`` doesn't expose a ``__dict__`` at all.
        """
        # En lugar de isinstance, verificamos que tenga los mismos atributos
        if not hasattr(other, "__dict__"):
            return "Invalid comparison. The object does not has the espected attributes."

        differences = {}
        for attr in vars(self):  # Obtener todos los atributos del objeto
            if getattr(self, attr) != getattr(other, attr):
                differences[attr] = (getattr(self, attr), getattr(other, attr))

        return differences if differences else "There are no changes."
