# pICkIT

Automated analysis of protein-ligand interactions.

`AnalyzeInteractions` parses Arpeggio or IChem interaction files into a
structured interaction matrix (`InteractionData`), then lets you filter,
sort, reshape, plot and export that matrix.

```python
from pickit import AnalyzeInteractions

analyzer = AnalyzeInteractions()
data = analyzer.analyze_files(directory="arpeggio_output", mode=analyzer.ARPEGGIO,
                               activity_file="activity.csv")
data = analyzer.sort_matrix(interaction_data=data)
analyzer.heatmap(interaction_data=data, title="Interactions", mode=analyzer.COUNT, save=True)
analyzer.save_interaction_data(interaction_data=data, filename="results.xlsx")
```

See the **API reference** in the sidebar for every public method, and
**Decisions** for the reasoning behind internal design choices made during
the [modularization](https://github.com/31ldts/pICkIT) of this project.

## Package structure

`AnalyzeInteractions` (in `pickit.analyzer`) is composed of five mixins,
each covering one concern:

| Mixin | Responsibility |
|---|---|
| `ValidationMixin` | shared input-validation helpers |
| `FilterMixin` | filtering, sorting, reshaping interaction matrices |
| `ExportMixin` | DataFrame / Excel export |
| `PlotMixin` | heatmaps, bar charts, pie charts |
| `IOMixin` | directory/config management, file parsing entry point |

File parsing itself lives in `pickit.parsers` (`arpeggio.py`, `ichem.py`,
`common.py`), called from `IOMixin.analyze_files`.

`pickit.analyze_interactions` is kept as a backward-compatibility shim:
`from pickit.analyze_interactions import AnalyzeInteractions` still works,
resolving to the exact same class as `from pickit import AnalyzeInteractions`.
