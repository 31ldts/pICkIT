import os

import pandas as pd
import pytest

from pickit.analyze_interactions import InvalidFileExtensionException, TypeMismatchException


class TestGetDataframe:
    def test_returns_dataframe_indexed_by_residue(self, analyzer, small_interaction_data):
        df = analyzer.get_dataframe(interaction_data=small_interaction_data)
        assert isinstance(df, pd.DataFrame)
        assert list(df.index) == ["HIS 41", "ASN 142", "GLU 166"]
        assert list(df.columns) == ["LIGA (7.5)", "LIGB (5.0)", "LIGC (0)"]

    def test_wrong_type_raises(self, analyzer):
        with pytest.raises(TypeMismatchException):
            analyzer.get_dataframe(interaction_data={"not": "an InteractionData"})


class TestGetInteractions:
    def test_returns_1_indexed_label_map(self, analyzer, small_interaction_data):
        result = analyzer.get_interactions(interaction_data=small_interaction_data)
        assert result[1] == "AMIDEAMIDE"
        assert result[13] == "hbond"
        assert len(result) == len(small_interaction_data.interactions)


class TestSaveInteractionData:
    def test_writes_xlsx_with_matrix_and_attributes_sheets(self, analyzer, small_interaction_data, tmp_path):
        analyzer.save_interaction_data(interaction_data=small_interaction_data, filename="out.xlsx")
        out_path = os.path.join(str(tmp_path), "out.xlsx")
        assert os.path.isfile(out_path)

        sheets = pd.read_excel(out_path, sheet_name=None)
        assert "Matrix" in sheets
        assert "Attributes" in sheets
        assert len(sheets["Attributes"]) == len(small_interaction_data.interactions)

    def test_non_xlsx_filename_raises(self, analyzer, small_interaction_data):
        with pytest.raises(InvalidFileExtensionException):
            analyzer.save_interaction_data(interaction_data=small_interaction_data, filename="out.csv")
