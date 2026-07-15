import pytest

from pickit.analyze_interactions import InvalidAxisException, MissedActivityException


class TestSortMatrix:
    def test_sorts_rows_by_interaction_count_descending(self, analyzer, small_interaction_data):
        out = analyzer.sort_matrix(interaction_data=small_interaction_data)
        row_labels = [row[0] for row in out.matrix[1:]]
        # ASN142 has 2 interactions in the fixture, HIS41 has 1, GLU166 has 0 -> ASN142 first
        assert row_labels[0] == "ASN 142"

    def test_thr_interactions_filters_low_count_rows(self, analyzer, small_interaction_data):
        out = analyzer.sort_matrix(interaction_data=small_interaction_data, thr_interactions=2)
        row_labels = [row[0] for row in out.matrix[1:]]
        assert row_labels == ["ASN 142"]

    def test_thr_activity_filters_by_ligand_activity(self, analyzer, small_interaction_data):
        out = analyzer.sort_matrix(interaction_data=small_interaction_data, thr_activity=6.0)
        header = out.matrix[0]
        # Only LIGA (7.5) clears the 6.0 threshold; LIGB (5.0) and LIGC (0) are dropped
        assert header == ["", "LIGA (7.5)"]

    def test_selected_items_keeps_top_n(self, analyzer, small_interaction_data):
        out = analyzer.sort_matrix(interaction_data=small_interaction_data, selected_items=1)
        row_labels = [row[0] for row in out.matrix[1:]]
        assert row_labels == ["ASN 142"]

    def test_conflicting_selection_criteria_raises(self, analyzer, small_interaction_data):
        with pytest.raises(ValueError):
            analyzer.sort_matrix(interaction_data=small_interaction_data, thr_interactions=1, selected_items=1)

    def test_invalid_axis_raises(self, analyzer, small_interaction_data):
        with pytest.raises(InvalidAxisException):
            analyzer.sort_matrix(interaction_data=small_interaction_data, axis="diagonal")

    def test_residue_chain_sorts_by_residue_number(self, analyzer, small_interaction_data):
        out = analyzer.sort_matrix(interaction_data=small_interaction_data, residue_chain=True)
        row_labels = [row[0] for row in out.matrix[1:]]
        assert row_labels == ["HIS 41", "ASN 142", "GLU 166"]


class TestRemoveEmptyAxis:
    def test_removes_empty_row_and_column(self, analyzer, small_interaction_data):
        out = analyzer.remove_empty_axis(interaction_data=small_interaction_data)
        row_labels = [row[0] for row in out.matrix[1:]]
        header = out.matrix[0]
        # GLU166 has no interactions at all -> row removed; LIGC (0) has no interactions -> column removed
        assert "GLU 166" not in row_labels
        assert "LIGC (0)" not in header

    def test_does_not_mutate_the_original(self, analyzer, small_interaction_data):
        original_rows = len(small_interaction_data.matrix)
        analyzer.remove_empty_axis(interaction_data=small_interaction_data)
        assert len(small_interaction_data.matrix) == original_rows


class TestTransposeMatrix:
    def test_swaps_rows_and_columns(self, analyzer, small_interaction_data):
        out = analyzer.transpose_matrix(interaction_data=small_interaction_data)
        assert len(out.matrix) == len(small_interaction_data.matrix[0])
        assert len(out.matrix[0]) == len(small_interaction_data.matrix)
        assert out.matrix[0][1] == "HIS 41"

    def test_double_transpose_is_identity(self, analyzer, small_interaction_data):
        out = analyzer.transpose_matrix(interaction_data=small_interaction_data)
        out = analyzer.transpose_matrix(interaction_data=out)
        assert out.matrix == small_interaction_data.matrix

    def test_too_small_matrix_raises(self, analyzer, small_interaction_data):
        small_interaction_data.matrix = [["only one cell"]]
        with pytest.raises(ValueError):
            analyzer.transpose_matrix(interaction_data=small_interaction_data)
