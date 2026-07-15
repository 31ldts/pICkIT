import os

import matplotlib
matplotlib.use("Agg")  # headless backend for CI / test environments

import pytest

from pickit.analyze_interactions import InvalidModeException


class TestHeatmap:
    def test_count_mode_saves_png(self, analyzer, small_interaction_data, tmp_path):
        analyzer.heatmap(interaction_data=small_interaction_data, title="count_heatmap", mode=analyzer.COUNT, save=True)
        assert os.path.isfile(os.path.join(str(tmp_path), "count_heatmap.png"))

    def test_max_mode_saves_png(self, analyzer, small_interaction_data, tmp_path):
        analyzer.heatmap(interaction_data=small_interaction_data, title="max_heatmap", mode=analyzer.MAXIMUM, save=True)
        assert os.path.isfile(os.path.join(str(tmp_path), "max_heatmap.png"))

    def test_invalid_mode_raises(self, analyzer, small_interaction_data):
        with pytest.raises(InvalidModeException):
            analyzer.heatmap(interaction_data=small_interaction_data, title="bad", mode="not_a_mode", save=True)


class TestBarChart:
    def test_saves_png(self, analyzer, small_interaction_data, tmp_path):
        analyzer.bar_chart(interaction_data=small_interaction_data, plot_name="bars", save=True)
        assert os.path.isfile(os.path.join(str(tmp_path), "bars.png"))

    def test_stacked_saves_png(self, analyzer, small_interaction_data, tmp_path):
        analyzer.bar_chart(interaction_data=small_interaction_data, plot_name="bars_stacked", stacked=True, save=True)
        assert os.path.isfile(os.path.join(str(tmp_path), "bars_stacked.png"))


class TestPieChart:
    def test_saves_png(self, analyzer, small_interaction_data, tmp_path):
        analyzer.pie_chart(interaction_data=small_interaction_data, plot_name="pie", axis=analyzer.ROWS, save=True)
        assert os.path.isfile(os.path.join(str(tmp_path), "pie.png"))

    def test_columns_axis_saves_png(self, analyzer, small_interaction_data, tmp_path):
        analyzer.pie_chart(interaction_data=small_interaction_data, plot_name="pie_cols", axis=analyzer.COLUMNS, save=True)
        assert os.path.isfile(os.path.join(str(tmp_path), "pie_cols.png"))
