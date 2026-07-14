"""Validation helpers, extracted as a mixin.

Extracted verbatim from the original monolithic ``analyze_interactions.py``
(Fase 2, paso 3 del plan de modularización). No logic changes were made
during this extraction — only the physical location of the code changed.

Implemented as a mixin (``ValidationMixin``) rather than free functions
because these methods are called as ``self._verify_dimensions(...)`` etc.
throughout the rest of the (still not fully split) class, and
``_verify_case`` reads ``self.LOWER`` / ``self.UPPER`` class constants
defined on ``AnalyzeInteractions``. Converting to free functions now would
mean touching every call site across the whole file in one shot, which
goes against the plan's "copiar-pegar primero, refactorizar después,
módulo a módulo" approach. Revisit this as a possible further
simplification once the whole class is mixin-composed (paso 8,
``analyzer.py``).
"""

from .exceptions import TypeMismatchException


class ValidationMixin:
    """Shared input-validation logic used by every other mixin."""

    def _check_variable_types(self, variables: list, expected_types: list, variable_names: list[str]):
        """
        Check if variables match their expected types.

        Args:
            variables (list): List of variables to check.
            expected_types (list): List of expected types (can include tuples).
            variable_names (list[str]): List of variable names for error messages.

        Raises:
            ValueError: If the lengths of input lists don't match.
            TypeMismatchException: If a variable type doesn't match the expected one.
        """
        if len(variables) != len(expected_types) or len(variables) != len(variable_names):
            raise ValueError(
                "The lists of variables, expected types, and variable names must all have the same length."
            )

        for i, variable in enumerate(variables):
            expected_type = expected_types[i]
            variable_name = variable_names[i]

            if not isinstance(variable, expected_type if isinstance(expected_type, tuple) else (expected_type,)):
                actual_type = type(variable)
                raise TypeMismatchException(
                    variable_name, expected_type if isinstance(expected_type, tuple) else (expected_type,), actual_type
                )

    def _get_residues_axis(self, matrix: list[list[str]]) -> str:
        """
        Determine whether residues' axis is in rows or columns.

        Args:
            matrix (list[list[str]]): Matrix with interaction data.

        Returns:
            str: 'columns' if residues are in columns, 'rows' otherwise.

        Raises:
            ValueError: If the axis cannot be determined.
        """
        self._verify_dimensions(matrix=matrix)

        # Check '(' occurrence in specific cells to determine axis
        count_0_1 = matrix[0][1].count("(")
        count_1_0 = matrix[1][0].count("(")

        if count_0_1 == count_1_0:  # This checks in case activity is used or not
            count_0_1 = matrix[0][1].count(", ")
            count_1_0 = matrix[1][0].count(", ")
            if count_0_1 == count_1_0:  # This checks in case PDB and ligand code are used or not
                count_0_1 = matrix[0][1].count(" ")
                count_1_0 = matrix[1][0].count(" ")
                if count_0_1 == count_1_0:
                    raise ValueError("Cannot determine the residues' axis.")
                elif count_0_1 == 1:
                    return "columns"
                elif count_1_0 == 1:
                    return "rows"
                else:
                    raise ValueError("Cannot determine the residues' axis.")
            elif count_0_1 == 1:
                return "rows"
            elif count_1_0 == 1:
                return "columns"
            else:
                raise ValueError("Cannot determine the residues' axis.")
        elif count_0_1 == 1:
            return "rows"
        elif count_1_0 == 1:
            return "columns"
        else:
            raise ValueError("Cannot determine the residues' axis.")

    def _verify_dimensions(self, matrix: list[list[str]], _pie: bool = False) -> None:
        """
        Verify matrix dimensions to ensure at least 2 rows and 2 columns.

        Args:
            matrix (list[list[str]]): Matrix to verify.

        Returns:
            None

        Raises:
            ValueError: If the matrix is too small or any row is too short.
        """
        if _pie:
            if len(matrix) < 1 or any(len(row) < 1 for row in matrix):
                raise ValueError("The input/output matrix is empty, it must have at least 1 row and 1 column.")
        elif len(matrix) < 2 or any(len(row) < 2 for row in matrix):
            raise ValueError("The input/output matrix is empty, it must have at least 2 rows and 2 columns.")

    def _verify_case(self, case: str) -> str:
        """
        Verifies the case parameter for heatmap visualization.

        Args:
            case (str): The case type to verify.

        Returns:
            str: The verified case type.

        Raises:
            ValueError: If the case is not recognized.
        """
        if case not in [self.LOWER, self.UPPER, None]:
            raise ValueError(f"Invalid case '{case}'. Expected 'lower', 'upper' or None.")
        return case
