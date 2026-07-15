import pytest

from pickit.analyze_interactions import TypeMismatchException


class TestFilterByInteraction:
    def test_keeps_only_requested_interaction_type(self, analyzer, small_interaction_data):
        out = analyzer.filter_by_interaction(interaction_data=small_interaction_data, interactions=[13])
        # hbond (13) survives in both cells that had it; hydrophobic (8) is stripped out
        assert out.matrix[1][1] == "13 |CA-O1(A)|"
        assert out.matrix[2][1] == "13 |CB-N1(A)|"
        assert out.matrix[2][2] == "13 |CB-N1(A)|"

    def test_does_not_mutate_the_original(self, analyzer, small_interaction_data):
        original = small_interaction_data.matrix[2][1]
        analyzer.filter_by_interaction(interaction_data=small_interaction_data, interactions=[13])
        assert small_interaction_data.matrix[2][1] == original

    def test_no_matching_interactions_raises_value_error(self, analyzer, small_interaction_data):
        with pytest.raises(ValueError):
            analyzer.filter_by_interaction(interaction_data=small_interaction_data, interactions=[1])

    def test_wrong_type_raises_type_mismatch(self, analyzer, small_interaction_data):
        with pytest.raises(TypeMismatchException):
            analyzer.filter_by_interaction(interaction_data=small_interaction_data, interactions="13")


class TestFilterBySubunit:
    def test_keeps_matching_subunit(self, analyzer, small_interaction_data):
        out = analyzer.filter_by_subunit(interaction_data=small_interaction_data, subunits=["A"])
        assert out.matrix == small_interaction_data.matrix

    def test_unknown_subunit_empties_cells(self, analyzer, small_interaction_data):
        out = analyzer.filter_by_subunit(interaction_data=small_interaction_data, subunits=["B"])
        # No interaction in the fixture belongs to subunit B, so every populated cell is emptied
        assert out.matrix[1][1] == "-"
        assert out.matrix[2][1] == "-"


class TestFilterByResidue:
    def test_filter_by_subpocket(self, analyzer, small_interaction_data):
        out = analyzer.filter_by_residue(
            interaction_data=small_interaction_data,
            subpocket_path="subpockets.csv",
            subpockets=["S1'"],
        )
        row_labels = [row[0] for row in out.matrix[1:]]
        # subpockets.csv defines S1' as containing HIS41 (among others); ASN142/GLU166 are not in S1'
        assert row_labels == ["HIS 41"]

    def test_filter_by_main_chain(self, analyzer, small_interaction_data):
        # CA is a main-chain atom; filtering by <main> should keep HIS41's CA-based interaction
        out = analyzer.filter_by_residue(interaction_data=small_interaction_data, chain="<main>")
        row_labels = [row[0] for row in out.matrix[1:]]
        assert "HIS 41" in row_labels

    def test_invalid_chain_raises(self, analyzer, small_interaction_data):
        from pickit.analyze_interactions import InvalidModeException
        with pytest.raises(InvalidModeException):
            analyzer.filter_by_residue(interaction_data=small_interaction_data, chain="<sideways>")
