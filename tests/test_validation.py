"""Tests for pickit._validation.ValidationMixin.

New tests introduced during Fase 2, step 3 (extraction of _validation.py).
These four methods are already exercised indirectly through almost every
other test in the suite (every public method calls into
`_check_variable_types` / `_verify_dimensions`), but had no dedicated,
isolated test file before. `_verify_case` in particular had no direct
baseline coverage at all (see Fase 0 notes on uncommon error branches).
"""

import pytest

from pickit._validation import ValidationMixin
from pickit.exceptions import TypeMismatchException


class _Host(ValidationMixin):
    """Minimal concrete host exposing the class constants _verify_case needs."""

    LOWER = "lower"
    UPPER = "upper"


@pytest.fixture
def host():
    return _Host()


class TestCheckVariableTypes:
    def test_passes_silently_when_types_match(self, host):
        host._check_variable_types(
            variables=[1, "a", [1, 2]],
            expected_types=[int, str, list],
            variable_names=["n", "s", "lst"],
        )  # no exception raised

    def test_accepts_tuple_of_expected_types(self, host):
        host._check_variable_types(variables=[None], expected_types=[(str, None.__class__)], variable_names=["x"])

    def test_raises_type_mismatch_with_details(self, host):
        with pytest.raises(TypeMismatchException) as exc_info:
            host._check_variable_types(variables=["ten"], expected_types=[int], variable_names=["plot_max_cols"])
        assert "plot_max_cols" in str(exc_info.value)

    def test_mismatched_list_lengths_raise_value_error(self, host):
        with pytest.raises(ValueError):
            host._check_variable_types(variables=[1, 2], expected_types=[int], variable_names=["a", "b"])


class TestVerifyDimensions:
    def test_accepts_2x2_or_larger(self, host):
        host._verify_dimensions(matrix=[["", "a"], ["b", "c"]])  # no exception

    def test_rejects_matrix_with_too_few_rows(self, host):
        with pytest.raises(ValueError):
            host._verify_dimensions(matrix=[["", "a"]])

    def test_rejects_row_that_is_too_short(self, host):
        with pytest.raises(ValueError):
            host._verify_dimensions(matrix=[["", "a"], ["b"]])

    def test_pie_mode_only_requires_1x1(self, host):
        host._verify_dimensions(matrix=[["a"]], _pie=True)  # no exception

    def test_pie_mode_rejects_empty_matrix(self, host):
        with pytest.raises(ValueError):
            host._verify_dimensions(matrix=[], _pie=True)


class TestGetResiduesAxis:
    def test_detects_rows_axis(self, host):
        # matrix[0][1] has 1 space (ligand code + activity), matrix[1][0] has 2 -> rows
        matrix = [["", "LIGA (7.5)"], ["HIS 41", ""]]
        assert host._get_residues_axis(matrix=matrix) == "rows"

    def test_detects_columns_axis(self, host):
        matrix = [["", "HIS 41"], ["LIGA (7.5)", ""]]
        assert host._get_residues_axis(matrix=matrix) == "columns"

    def test_undeterminable_axis_raises(self, host):
        matrix = [["", "AAA"], ["BBB", ""]]
        with pytest.raises(ValueError):
            host._get_residues_axis(matrix=matrix)


class TestVerifyCase:
    def test_accepts_lower(self, host):
        assert host._verify_case("lower") == "lower"

    def test_accepts_upper(self, host):
        assert host._verify_case("upper") == "upper"

    def test_accepts_none(self, host):
        assert host._verify_case(None) is None

    def test_rejects_anything_else(self, host):
        with pytest.raises(ValueError):
            host._verify_case("mixed")


def test_analyze_interactions_inherits_validation_mixin():
    from pickit.analyze_interactions import AnalyzeInteractions

    assert issubclass(AnalyzeInteractions, ValidationMixin)
