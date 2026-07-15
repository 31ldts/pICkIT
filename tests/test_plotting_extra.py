"""Supplementary plot_mixin.py tests, added during Fase 3.

The Fase 0 baseline explicitly flagged these branches as uncovered
("Cobertura no ejercida principalmente en... buena parte del cuerpo
interno de heatmap/bar_chart/pie_chart que solo se ejerce con datasets
grandes... o con subpocket_path/export_xlsx" — see plan-accion-modularizacion-pickit.md).
These tests close that gap using the real project fixtures
(tests/data/subpockets.csv, tests/data/mini_activity.csv), without
touching the frozen Fase 0 test files.
"""

import os

import matplotlib
import pandas as pd

matplotlib.use("Agg")  # headless backend for CI / test environments

import pytest

from pickit.analyze_interactions import MissedActivityException


class TestHeatmapActivityModes:
    def test_max_mode_with_real_activity_saves_png(self, analyzer, arpeggio_dir, activity_csv, tmp_path):
        data = analyzer.analyze_files(directory=arpeggio_dir, mode=analyzer.ARPEGGIO, activity_file=activity_csv)
        analyzer.heatmap(interaction_data=data, title="max_activity_heatmap", mode=analyzer.MAXIMUM, save=True)
        assert os.path.isfile(os.path.join(str(tmp_path), "max_activity_heatmap.png"))

    def test_mean_mode_with_real_activity_saves_png(self, analyzer, arpeggio_dir, activity_csv, tmp_path):
        data = analyzer.analyze_files(directory=arpeggio_dir, mode=analyzer.ARPEGGIO, activity_file=activity_csv)
        analyzer.heatmap(interaction_data=data, title="mean_activity_heatmap", mode=analyzer.MEAN, save=True)
        assert os.path.isfile(os.path.join(str(tmp_path), "mean_activity_heatmap.png"))

    def test_min_max_mean_without_activity_file_raises(self, analyzer, arpeggio_dir):
        # No activity_file passed -> ligand columns have no "(value)" suffix,
        # so min/max/mean can't parse an activity and must raise.
        data = analyzer.analyze_files(directory=arpeggio_dir, mode=analyzer.ARPEGGIO)
        with pytest.raises(MissedActivityException):
            analyzer.heatmap(interaction_data=data, title="no_activity", mode=analyzer.MAXIMUM, save=True)


class TestHeatmapSubpocketColoring:
    def test_heatmap_with_subpocket_path_saves_png(self, analyzer, small_interaction_data, tmp_path):
        analyzer.heatmap(
            interaction_data=small_interaction_data,
            title="subpocket_heatmap",
            mode=analyzer.COUNT,
            subpocket_path="subpockets.csv",
            save=True,
        )
        assert os.path.isfile(os.path.join(str(tmp_path), "subpocket_heatmap.png"))

    def test_heatmap_with_custom_subpocket_colors_saves_png(self, analyzer, small_interaction_data, tmp_path):
        analyzer.heatmap(
            interaction_data=small_interaction_data,
            title="subpocket_custom_colors",
            mode=analyzer.COUNT,
            subpocket_path="subpockets.csv",
            subpocket_colors=["#8C2CE2", "#DFBA52", "#00AA88", "#FF5533"],
            save=True,
        )
        assert os.path.isfile(os.path.join(str(tmp_path), "subpocket_custom_colors.png"))

    def test_heatmap_split_by_atom_saves_png(self, analyzer, small_interaction_data, tmp_path):
        analyzer.heatmap(
            interaction_data=small_interaction_data,
            title="split_by_atom_heatmap",
            mode=analyzer.COUNT,
            split_by_atom=True,
            save=True,
        )
        assert os.path.isfile(os.path.join(str(tmp_path), "split_by_atom_heatmap.png"))

    def test_heatmap_remove_empty_drops_all_nan_rows(self, analyzer, small_interaction_data, tmp_path):
        analyzer.heatmap(
            interaction_data=small_interaction_data,
            title="remove_empty_heatmap",
            mode=analyzer.COUNT,
            remove_empty=True,
            save=True,
        )
        assert os.path.isfile(os.path.join(str(tmp_path), "remove_empty_heatmap.png"))


class TestExportXlsxBranches:
    def test_bar_chart_export_xlsx_not_stacked(self, analyzer, small_interaction_data, tmp_path):
        analyzer.bar_chart(interaction_data=small_interaction_data, plot_name="bars_xlsx", save=True, export_xlsx=True)
        xlsx_path = os.path.join(str(tmp_path), "bars_xlsx.xlsx")
        assert os.path.isfile(xlsx_path)
        df = pd.read_excel(xlsx_path)
        assert "Element" in df.columns

    def test_bar_chart_export_xlsx_stacked(self, analyzer, small_interaction_data, tmp_path):
        analyzer.bar_chart(
            interaction_data=small_interaction_data,
            plot_name="bars_xlsx_stacked",
            stacked=True,
            save=True,
            export_xlsx=True,
        )
        xlsx_path = os.path.join(str(tmp_path), "bars_xlsx_stacked.xlsx")
        assert os.path.isfile(xlsx_path)

    def test_pie_chart_export_xlsx(self, analyzer, small_interaction_data, tmp_path):
        analyzer.pie_chart(
            interaction_data=small_interaction_data,
            plot_name="pie_xlsx",
            axis=analyzer.ROWS,
            save=True,
            export_xlsx=True,
        )
        xlsx_path = os.path.join(str(tmp_path), "pie_xlsx.xlsx")
        assert os.path.isfile(xlsx_path)
        df = pd.read_excel(xlsx_path)
        assert "Interaction" in df.columns and "Percent" in df.columns


class TestBarChartSubpocketColoring:
    def test_bar_chart_with_subpocket_path_saves_png(self, analyzer, small_interaction_data, tmp_path):
        analyzer.bar_chart(
            interaction_data=small_interaction_data,
            plot_name="bars_subpocket",
            subpocket_path="subpockets.csv",
            save=True,
        )
        assert os.path.isfile(os.path.join(str(tmp_path), "bars_subpocket.png"))


class TestCaseFormatting:
    def test_heatmap_upper_case_labels(self, analyzer, small_interaction_data, tmp_path):
        analyzer.heatmap(
            interaction_data=small_interaction_data,
            title="upper_case",
            mode=analyzer.COUNT,
            case=analyzer.UPPER,
            save=True,
        )
        assert os.path.isfile(os.path.join(str(tmp_path), "upper_case.png"))

    def test_bar_chart_lower_case_labels(self, analyzer, small_interaction_data, tmp_path):
        analyzer.bar_chart(
            interaction_data=small_interaction_data,
            plot_name="bars_lower",
            stacked=True,
            case=analyzer.LOWER,
            save=True,
        )
        assert os.path.isfile(os.path.join(str(tmp_path), "bars_lower.png"))
