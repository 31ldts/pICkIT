"""
Protein-ligand interaction analysis for the SARS-CoV-2 Mpro dataset (386 complexes).

Real usage example of pICkIT: parses the Arpeggio dataset annotated with activity (pIC50)
and answers 4 independent research questions, each individually toggleable
in QUESTIONS_TO_RUN below.
"""

from pickit.analyze_interactions import AnalyzeInteractions

# --------------------------------------------------------------------------- #
# Configuration
# --------------------------------------------------------------------------- #

INPUT_DIR = "input"
OUTPUT_DIR = "output"
ARPEGGIO_DIR = "386_arpeggio"
ACTIVITY_FILE = "386_Mpro_nc_February_2025.csv"
TEMPLATE_FILE = "template.json"
SUBPOCKETS_FILE = "subpockets.csv"
RAW_MATRIX_FILE = "prueba.xlsx"

ACTIVITY_THRESHOLD = 7.4         # minimum pIC50 to consider a compound "active"
TOP_RESIDUES_THRESHOLD = 200     # minimum number of interactions to consider a residue "top"

# Which questions to run in this pass
QUESTIONS_TO_RUN = {
    "q1_most_frequent_interactions": True,
    "q2_activity_vs_interaction_type": True,
    "q3_interaction_count_vs_activity": True,
    "q4_subpocket_occupancy": True,
}


# --------------------------------------------------------------------------- #
# Setup
# --------------------------------------------------------------------------- #

def build_analyzer() -> AnalyzeInteractions:
    """Creates the analyzer and points it to the input/output directories."""
    analyzer = AnalyzeInteractions()
    analyzer.change_directory(OUTPUT_DIR, mode=analyzer.OUTPUT)
    analyzer.change_directory(INPUT_DIR, mode=analyzer.INPUT)
    analyzer.set_config(heat_max_cols=10)
    return analyzer


def load_data(analyzer: AnalyzeInteractions):
    """Parses the Arpeggio dataset annotated with activity and saves the raw matrix."""
    return analyzer.analyze_files(
        directory=ARPEGGIO_DIR,
        mode=analyzer.ARPEGGIO,
        activity_file=ACTIVITY_FILE,
        template_file=TEMPLATE_FILE,
        save=RAW_MATRIX_FILE,
    )


# --------------------------------------------------------------------------- #
# Q1 — Which interactions and residues are most frequent?
# --------------------------------------------------------------------------- #

def q1_most_frequent_interactions(analyzer: AnalyzeInteractions, data) -> None:
    analyzer.pie_chart(
        interaction_data=data,
        plot_name="Q1 Baseline interaction type distribution",
        title="Baseline interaction type distribution",
        axis=analyzer.ROWS,
        save=True,
    )
    analyzer.bar_chart(
        interaction_data=data,
        plot_name="Q1 Baseline residue-interaction profile",
        title="Baseline residue-interaction profile",
        axis=analyzer.ROWS,
        stacked=True,
        save=True,
        type_count=False,
    )

    top_residues = analyzer.sort_matrix(
        interaction_data=data, thr_interactions=TOP_RESIDUES_THRESHOLD, axis=analyzer.ROWS
    )
    analyzer.pie_chart(
        interaction_data=top_residues,
        plot_name="Q1 Top residues interaction type distribution",
        title="Top residues interaction type distribution",
        axis=analyzer.ROWS,
        save=True,
    )
    analyzer.bar_chart(
        interaction_data=top_residues,
        plot_name="Q1 Top residues interaction profile",
        title="Top residues interaction profile",
        axis=analyzer.ROWS,
        stacked=True,
        save=True,
        type_count=False,
    )


# --------------------------------------------------------------------------- #
# Q2 — Which interactions predominate in the most active compounds (high pIC50)?
# --------------------------------------------------------------------------- #

def q2_activity_vs_interaction_type(analyzer: AnalyzeInteractions, data) -> None:
    heatmap_modes = (analyzer.MAXIMUM, analyzer.MEAN, analyzer.COUNT)
    for mode in heatmap_modes:
        top = analyzer.sort_matrix(interaction_data=data, thr_activity=ACTIVITY_THRESHOLD, axis=analyzer.ROWS)
        top = analyzer.remove_empty_axis(interaction_data=top)
        analyzer.heatmap(interaction_data=top, title="", mode=mode, save=False)


# --------------------------------------------------------------------------- #
# Q3 — Is there a correlation between the number of interactions (total and by type) and activity?
# --------------------------------------------------------------------------- #

def q3_interaction_count_vs_activity(analyzer: AnalyzeInteractions, data) -> None:
    analyzer.pie_chart(
        interaction_data=data,
        plot_name="Q3 Baseline interaction type distribution",
        title="Baseline interaction type distribution",
        axis=analyzer.ROWS,
        save=False,
    )
    active_compounds = analyzer.sort_matrix(interaction_data=data, thr_activity=ACTIVITY_THRESHOLD)
    analyzer.pie_chart(
        interaction_data=active_compounds,
        plot_name=f"Q3-InteMasFrecPie({ACTIVITY_THRESHOLD})",
        axis=analyzer.ROWS,
        save=False,
    )


# --------------------------------------------------------------------------- #
# Q4 — Which Mpro subsites does each compound occupy?
# --------------------------------------------------------------------------- #

def q4_subpocket_occupancy(analyzer: AnalyzeInteractions, data) -> None:
    subpockets = ["S1'", "S1", "S2", "S4"]
    filtered_by_subpocket = {}

    for subpocket in subpockets:
        subset = analyzer.filter_by_residue(
            interaction_data=data, subpocket_path=SUBPOCKETS_FILE, subpockets=[subpocket]
        )
        safe_name = subpocket.replace("'", "")
        subset = analyzer.remove_empty_axis(interaction_data=subset, save=f"Q4 {safe_name}.xlsx")
        filtered_by_subpocket[subpocket] = subset

    subpocket_to_detail = filtered_by_subpocket["S4"]

    analyzer.heatmap(interaction_data=subpocket_to_detail, title="A", mode=analyzer.MAXIMUM, save=False)
    analyzer.heatmap(interaction_data=subpocket_to_detail, title="B", mode=analyzer.MEAN, save=False)
    analyzer.heatmap(interaction_data=subpocket_to_detail, title="C", mode=analyzer.COUNT, save=False)

    analyzer.pie_chart(
        interaction_data=subpocket_to_detail,
        plot_name="Q4 Baseline interaction type distribution (S1')",
        title="Baseline interaction type distribution (S1')",
        axis=analyzer.ROWS,
        save=True,
    )


# --------------------------------------------------------------------------- #
# Entry point
# --------------------------------------------------------------------------- #

def main() -> None:
    analyzer = build_analyzer()
    data = load_data(analyzer)

    print(analyzer.get_interactions(data))

    if QUESTIONS_TO_RUN["q1_most_frequent_interactions"]:
        q1_most_frequent_interactions(analyzer, data)
    if QUESTIONS_TO_RUN["q2_activity_vs_interaction_type"]:
        q2_activity_vs_interaction_type(analyzer, data)
    if QUESTIONS_TO_RUN["q3_interaction_count_vs_activity"]:
        q3_interaction_count_vs_activity(analyzer, data)
    if QUESTIONS_TO_RUN["q4_subpocket_occupancy"]:
        q4_subpocket_occupancy(analyzer, data)


if __name__ == "__main__":
    main()
