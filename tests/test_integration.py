"""End-to-end integration test for AnalyzeInteractions.

Fase 3 checklist item: a single test exercising the full pipeline
(parse -> filter/sort -> plot -> export) exactly as a user following the
README/quick start would, to verify the 5-mixin composition
(ValidationMixin, FilterMixin, ExportMixin, PlotMixin, IOMixin) works
together as a whole, not just as individually-tested pieces.

Uses the same tiny synthetic Arpeggio fixtures as the rest of the suite
(tests/data/arpeggio_mini, tests/data/mini_activity.csv), never the real
365 MB datasets, per the Fase 0 rule.
"""

import os

import matplotlib

matplotlib.use("Agg")  # headless backend for CI / test environments

import pandas as pd

from pickit import AnalyzeInteractions


def test_full_pipeline_parse_filter_sort_plot_export(tmp_path):
    analyzer = AnalyzeInteractions()
    analyzer.input_directory = os.path.join(os.path.dirname(__file__), "data")
    analyzer.saving_directory = str(tmp_path)

    # 1. Parse (io_mixin.analyze_files)
    data = analyzer.analyze_files(
        directory="arpeggio_mini",
        mode=analyzer.ARPEGGIO,
        activity_file="mini_activity.csv",
    )
    assert len(data.matrix) > 1
    assert len(data.matrix[0]) > 1

    # 2. Filter + sort + reshape (filter_mixin)
    sorted_data = analyzer.sort_matrix(interaction_data=data)
    filtered_data = analyzer.filter_by_subunit(interaction_data=sorted_data, subunits=["A"])
    cleaned_data = analyzer.remove_empty_axis(interaction_data=filtered_data)
    assert len(cleaned_data.matrix) >= 2
    assert len(cleaned_data.matrix[0]) >= 2

    # 3. Plot (plot_mixin) — just confirm it renders and saves without raising
    analyzer.heatmap(interaction_data=cleaned_data, title="integration_heatmap", mode=analyzer.COUNT, save=True)
    assert os.path.isfile(os.path.join(str(tmp_path), "integration_heatmap.png"))

    analyzer.bar_chart(interaction_data=cleaned_data, plot_name="integration_bars", save=True)
    assert os.path.isfile(os.path.join(str(tmp_path), "integration_bars.png"))

    # 4. Export (export_mixin)
    df = analyzer.get_dataframe(interaction_data=cleaned_data)
    assert isinstance(df, pd.DataFrame)
    assert not df.empty

    analyzer.save_interaction_data(interaction_data=cleaned_data, filename="integration_output.xlsx")
    out_path = os.path.join(str(tmp_path), "integration_output.xlsx")
    assert os.path.isfile(out_path)

    sheets = pd.read_excel(out_path, sheet_name=None)
    assert "Matrix" in sheets
    assert "Attributes" in sheets


def test_same_analyzer_instance_reusable_across_multiple_files(tmp_path):
    """A single AnalyzeInteractions instance can run the pipeline twice
    without leaking state between calls (config reset via set_config,
    fresh InteractionData objects each time)."""
    analyzer = AnalyzeInteractions()
    analyzer.input_directory = os.path.join(os.path.dirname(__file__), "data")
    analyzer.saving_directory = str(tmp_path)

    first = analyzer.analyze_files(directory="arpeggio_mini", mode=analyzer.ARPEGGIO)
    second = analyzer.analyze_files(directory="arpeggio_mini", mode=analyzer.ARPEGGIO)

    assert first.matrix == second.matrix
    assert first is not second
