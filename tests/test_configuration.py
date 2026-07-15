import os

import pytest

from pickit.analyze_interactions import (
    AnalyzeInteractions,
    InvalidModeException,
    InvalidColorException,
    TypeMismatchException,
)


class TestChangeDirectory:
    def test_happy_path_sets_input_directory(self, tmp_path):
        a = AnalyzeInteractions()
        sub = tmp_path / "input"
        sub.mkdir()
        cwd = os.getcwd()
        try:
            os.chdir(tmp_path)
            a.change_directory("input", mode=a.INPUT)
            assert a.input_directory == os.path.join(str(tmp_path), "input")
        finally:
            os.chdir(cwd)

    def test_happy_path_sets_output_directory(self, tmp_path):
        a = AnalyzeInteractions()
        sub = tmp_path / "output"
        sub.mkdir()
        cwd = os.getcwd()
        try:
            os.chdir(tmp_path)
            a.change_directory("output", mode=a.OUTPUT)
            assert a.saving_directory == os.path.join(str(tmp_path), "output")
        finally:
            os.chdir(cwd)

    def test_nonexistent_directory_raises_value_error(self, tmp_path):
        a = AnalyzeInteractions()
        cwd = os.getcwd()
        try:
            os.chdir(tmp_path)
            with pytest.raises(ValueError):
                a.change_directory("does_not_exist", mode=a.INPUT)
        finally:
            os.chdir(cwd)

    def test_invalid_mode_raises_invalid_mode_exception(self, tmp_path):
        a = AnalyzeInteractions()
        cwd = os.getcwd()
        try:
            os.chdir(tmp_path)
            with pytest.raises(InvalidModeException):
                a.change_directory(str(tmp_path), mode="sideways")
        finally:
            os.chdir(cwd)

    def test_wrong_type_raises_type_mismatch(self):
        a = AnalyzeInteractions()
        with pytest.raises(TypeMismatchException):
            a.change_directory(123, mode=a.INPUT)


class TestSetConfig:
    def test_default_config(self):
        a = AnalyzeInteractions()
        assert a.plot_max_cols == 80
        assert a.heat_max_cols == 30
        assert a.heat_colors == "RdYlGn"

    def test_arpeggio_mode_sets_labels_and_colors(self):
        a = AnalyzeInteractions()
        a.set_config(mode="arpeggio")
        assert "hbond" in a.interaction_labels
        assert len(a.plot_colors) > 0

    def test_ichem_mode_resets_to_default_labels(self):
        a = AnalyzeInteractions()
        a.set_config(mode="arpeggio")
        a.set_config(mode="ichem")
        assert "Hydrophobic" in a.interaction_labels

    def test_invalid_mode_raises(self):
        a = AnalyzeInteractions()
        with pytest.raises(InvalidModeException):
            a.set_config(mode="not_a_mode")

    def test_custom_plot_max_cols(self):
        a = AnalyzeInteractions()
        a.set_config(plot_max_cols=5)
        assert a.plot_max_cols == 5

    def test_heat_max_cols_and_colors(self):
        a = AnalyzeInteractions()
        a.set_config(heat_max_cols=10, heat_colors="viridis")
        assert a.heat_max_cols == 10
        assert a.heat_colors == "viridis"

    def test_invalid_hex_color_raises(self):
        a = AnalyzeInteractions()
        with pytest.raises(InvalidColorException):
            a.set_config(plot_colors=["not-a-color"])

    def test_valid_hex_colors_accepted(self):
        a = AnalyzeInteractions()
        a.set_config(plot_colors=["#ffffff", "#000000"])
        assert a.plot_colors == ["#ffffff", "#000000"]

    def test_reset_restores_defaults(self):
        a = AnalyzeInteractions()
        a.set_config(plot_max_cols=5, heat_max_cols=1)
        a.set_config(reset=True)
        assert a.plot_max_cols == 80
        assert a.heat_max_cols == 30

    def test_wrong_type_raises_type_mismatch(self):
        a = AnalyzeInteractions()
        with pytest.raises(TypeMismatchException):
            a.set_config(plot_max_cols="ten")
