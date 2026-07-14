"""Plotting operations: heatmaps, bar charts and pie charts.

Extracted verbatim from the original monolithic ``analyze_interactions.py``
(Fase 2, paso 6 del plan de modularización). No logic changes were made
during this extraction — only the physical location of the code changed.
This is the largest and highest-risk module of the whole modularization
after the parsers, so it's extracted as a single block rather than
piecemeal, exactly as the plan orders it.
"""

import copy
import operator
import os

import matplotlib.colors as mcolors
import mplcursors
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import colormaps
from matplotlib import pyplot as plt
from matplotlib.patches import Patch
from matplotlib.ticker import MaxNLocator

from .constants import DIFF_DELIM, GROUP_DELIM, HEATMAP_MODES, is_not_empty_or_dash
from .exceptions import HeatmapActivityException, InvalidModeException, MissedActivityException
from .models import InteractionData


class PlotMixin:
    """Heatmap, bar-chart and pie-chart generation.

    Depends on ``ValidationMixin`` (``self._verify_case``), ``FilterMixin``
    (``self.transpose_matrix``, ``self._stack_reactives``), ``ExportMixin``
    (``self._export_bar_data_to_excel``, ``self._export_pie_data_to_excel``),
    and on configuration state (``self.plot_colors``, ``self.plot_max_cols``,
    ``self.heat_max_cols``, ``self.heat_colors``, ``self.interaction_labels``,
    ``self.saving_directory``, ``self.input_directory``, ``self.set_config``).
    """

    def _plot_init(self, colors, matrix, axis, type_count):
        if colors is None:
            colors = self.plot_colors

        # Ensure the number of colors matches the number of interaction labels
        if len(colors) < len(self.interaction_labels):
            raise ValueError(
                f"Not enough colors provided. Expected at least {len(self.interaction_labels)} colors, "
                f"but got {len(colors)}."
            )

        # Calculate stacked data if necessary
        data, indices = self._stack_reactives(matrix=matrix, axis=axis, type_count=type_count)
        transposed_data = self.transpose_matrix(data, _pie=True)
        return colors, data, indices, transposed_data

    def _plot_end(self, save, plt, fig, plot_name):
        # Show or save the plot
        if not save:
            plt.show()
        else:
            plt.savefig(os.path.join(self.saving_directory, plot_name + ".png"))
            plt.close(fig)  # Close the figure after saving to avoid display overlap

    def _load_subsite_color_maps(self, csv_path: str, subsite_colors: list[str] = None):
        """
        Load subsite definitions from a CSV file and generate two color maps:

        1) residue_color_map : maps "GLU 166" → "#RRGGBB"
        2) subsite_color_map : maps "S1" → "#RRGGBB"

        Rules:
        - CSV format: Subsite, "GLU166, SER144<main>, VAL42<side>, ..."
        - Removes <main> or <side>
        - Converts "GLU166" → "GLU 166"
        - If fewer than 2 subsites exist → return empty dicts
        - If custom subsite_colors list is provided:
            * If too short → colors repeat cyclically
            * If too long → excess ignored
        """
        import csv

        # --- 1. Read CSV into a dict: subsite → residues ---
        csv_path = os.path.join(self.input_directory, csv_path)
        subsite_to_residues = {}

        try:
            with open(csv_path, mode="r", encoding="utf-8") as f:
                reader = csv.reader(f)
                for row in reader:
                    if not row or len(row) < 2:
                        continue

                    subsite = row[0].strip()
                    residues_raw = row[1].split(",")

                    cleaned_res = []
                    for item in residues_raw:
                        res = item.strip()
                        side_chain = None  # Variable para almacenar si es main o side

                        # Detectar y extraer <main> o <side>
                        if "<" in res:
                            # Extraer lo que está dentro de <>
                            if "<main>" in res:
                                side_chain = "main"
                                res = res.replace("<main>", "").strip()
                            elif "<side>" in res:
                                side_chain = "side"
                                res = res.replace("<side>", "").strip()
                            else:
                                # Si hay otro <...> que no sea main/side, solo quitarlo
                                res = res.split("<")[0].strip()

                        # Convertir "GLU166" → "GLU 166" (si no hay espacio ya)
                        # Primero verificamos el formato sin espacios
                        if len(res) > 3 and res[:3].isalpha() and res[3:].split(" ")[0].replace(" ", "").isdigit():
                            # Asegurar un máximo de un espacio en la parte numérica
                            parts = res[:3] + " " + res[3:]
                            # Normalizar: eliminar espacios múltiples y dejar solo uno
                            parts = " ".join(parts.split())
                            res = parts

                        # Aplicar modificaciones según side_chain
                        if side_chain == "main":
                            # Añadir los átomos de la cadena principal
                            # Asumiendo átomos estándar: N, CA, C, O
                            for atom in ["N", "CA", "C", "O"]:
                                cleaned_res.append(f"{res} {atom}")
                            continue
                        elif side_chain == "side":
                            # Añadir un + al final
                            res = f"{res} +"

                        cleaned_res.append(res)

                    subsite_to_residues[subsite] = cleaned_res

        except Exception as e:
            print("[WARNING] Failed to load subsite CSV:", e)
            return {}, {}

        # --- 2. Fewer than 2 subsites → no coloring ---
        if len(subsite_to_residues) < 2:
            return {}, {}

        # --- 3. Assign colors to subsites ---
        n_subsites = len(subsite_to_residues)

        if subsite_colors is None:
            cmap = colormaps.get_cmap("Dark2").resampled(n_subsites)
            final_colors = [mcolors.to_hex(cmap(i)) for i in range(n_subsites)]
        elif isinstance(subsite_colors, str):
            try:
                cmap = colormaps.get_cmap(subsite_colors).resampled(n_subsites)
                final_colors = [mcolors.to_hex(cmap(i)) for i in range(n_subsites)]
            except KeyError:
                print(f"[WARNING] Unknown colormap '{subsite_colors}'. Falling back to 'Dark2'.")
                cmap = colormaps.get_cmap("Dark2").resampled(n_subsites)
                final_colors = [mcolors.to_hex(cmap(i)) for i in range(n_subsites)]
        elif isinstance(subsite_colors, list):
            # Ensure color list is long enough
            if len(subsite_colors) < n_subsites:
                needed = n_subsites - len(subsite_colors)
                subsite_colors = subsite_colors + subsite_colors[:needed]
            final_colors = subsite_colors[:n_subsites]
        else:
            print("[WARNING] subsite_colors must be None, a list, or a string. Using default palette.")
            cmap = colormaps.get_cmap("Dark2").resampled(n_subsites)
            final_colors = [mcolors.to_hex(cmap(i)) for i in range(n_subsites)]

        # Build subsite → color dict
        subsite_color_map = {subsite: final_colors[i] for i, subsite in enumerate(subsite_to_residues.keys())}

        # --- 4. Build residue → color dict ---
        residue_color_map = {}
        for subsite, residues in subsite_to_residues.items():
            col = subsite_color_map[subsite]
            for r in residues:
                residue_color_map[r] = col

        return residue_color_map, subsite_color_map

    def _get_contrasting_text_color(self, bg_color):
        """
        Return 'black' or 'white' depending on which one contrasts better
        with the given background color.

        This works for:
        - CSS color names ('red', 'navy', ...)
        - Hex colors ('#RRGGBB')
        - RGB/RGBA tuples (0-1 float values)

        Logic:
        - Convert background color to RGB
        - Compute luminance using WCAG formula
        - Return white text for dark backgrounds, black text for light backgrounds
        """

        try:
            rgb = mcolors.to_rgb(bg_color)
        except ValueError:
            # Fallback: unknown color → assume white background
            return "black"

        # Compute luminance (per W3C accessibility standard)
        luminance = 0.2126 * rgb[0] + 0.7152 * rgb[1] + 0.0722 * rgb[2]

        # Threshold: 0.6 gives good contrast in practice
        return "black" if luminance > 0.6 else "white"

    def heatmap(
        self,
        interaction_data: InteractionData,
        title: str,
        mode: str,
        x_label: str = "",
        y_label: str = "",
        min_v: int = None,
        max_v: int = None,
        case: str = None,
        subpocket_path: str = None,
        subpocket_colors: list[str] = None,
        remove_empty: bool = False,
        split_by_atom: bool = False,
        save: bool = False,
    ):
        """
        Generates a heatmap based on interaction data using different processing modes.

        This method processes the interaction matrix, computes interaction statistics
        according to the specified mode, and generates a heatmap visualization. If the
        dataset is large, the visualization is split into multiple heatmaps.

        Supported processing modes:
            - 'min': Displays the minimum interaction values.
            - 'max': Displays the maximum interaction values.
            - 'mean': Computes and visualizes the average interaction values.
            - 'count': Counts the occurrences of interactions.
            - 'percent': Displays the percentage of interactions.

        Args:
            interaction_data (InteractionData): The object containing the interaction matrix.
            title (str): Title of the heatmap.
            mode (str): Processing mode ('min', 'max', 'mean', 'count', or 'percent').
            x_label (str, optional): Label for the x-axis. Defaults to an empty string.
            y_label (str, optional): Label for the y-axis. Defaults to an empty string.
            min_v (int, optional): Minimum value for the heatmap color scale. Defaults to None (auto-scaling).
            max_v (int, optional): Maximum value for the heatmap color scale. Defaults to None (auto-scaling).
            case (str, optional): Case style for the heatmap labels. Can be 'upper', 'lower', or None. Defaults to None.
            save (bool, optional): If True, saves the heatmap instead of displaying it. Defaults to False.
            subpocket_path (str, optional): Path to a CSV defining subpockets and their residues.
               Used to color residue labels by subsite. If the file defines fewer than two
               subsites, coloring is disabled.
            subpocket_colors (list[str] or str, optional): Custom colors for subpockets
               (e.g., ['#8C2CE2','#DFBA52'], colors are used in order and repeated if needed),
               or a colormap name (e.g., "Dark2", "tab10". it is interpreted as a Matplotlib colormap name).

        Returns:
            None: The function either displays the heatmap(s) or saves them to files.

        Raises:
            InvalidModeException: If an unsupported processing mode is provided.
            HeatmapActivityException: If the matrix contains invalid or negative activity values.
        """
        residue_color_map = None
        subsite_color_map = None

        if subpocket_path:
            residue_color_map, subsite_color_map = self._load_subsite_color_maps(subpocket_path, subpocket_colors)

        def validate_and_prepare_matrix(matrix: list[list[str]], mode: str):
            """
            Validates the input matrix and prepares it for processing.

            This method checks that all ligand activity values in the matrix are non-negative and transposes the matrix
            if the residue axis is configured to be in columns.

            Args:
                matrix (list[list[str]]): The matrix to validate and prepare.

            Returns:
                list[list[str]]: The validated matrix.

            Raises:
                HeatmapActivityException: If any ligand activity value is negative.
            """
            if self._get_residues_axis == "columns":
                matrix = self.transpose_matrix(matrix)
            if mode in ["min", "max", "mean"]:
                for ligand in matrix[0][1:]:
                    activity = ligand.split("(")[-1].replace(")", "")
                    try:
                        if float(activity) < 0.0:
                            raise HeatmapActivityException
                    except ValueError:
                        raise MissedActivityException("The matrix does not contain activity values.")
            return matrix

        def process_matrix(
            matrix: list[list[str]], mode: str, nan_template: list[float], split_by_atom: bool = False
        ) -> dict:
            """
            Procesa la matriz según el modo indicado, creando dinámicamente las entradas
            del diccionario de datos. Si split_by_atom es True, desglosa cada interacción
            por los átomos concretos del residuo que participan.

            Args:
                matrix (list[list[str]]): Matriz transpuesta (filas = residuos, columnas = ligandos).
                mode (str): Modo de procesamiento ('min','max','mean','count','percent').
                nan_template (list[float]): Plantilla de NaN para inicializar nuevos residuos.
                split_by_atom (bool): Si es True, crea una entrada por cada átomo del residuo
                                    que interviene en la interacción.

            Returns:
                dict: Diccionario {residuo (o residuo átomo): lista de valores por tipo de interacción}.
            """
            data = {}
            op = {"min": operator.lt, "max": operator.gt}.get(mode, None)

            # matrix[0] contiene la cabecera (residuos), el resto son filas de residuos
            for line in matrix[1:]:
                # La actividad se extrae de la primera columna de cada fila (nombre del ligando)
                if mode in ("min", "max", "mean"):
                    activity = float(line[0].split("(")[-1].replace(")", ""))

                for index in range(1, len(line)):
                    residue_base = matrix[0][index].split("-")[0]  # nombre base del residuo
                    cell = line[index]
                    sections = cell.split(DIFF_DELIM)

                    for section in sections:
                        # Ignorar celdas vacías o guiones
                        if not is_not_empty_or_dash(section.split(" ")[0]):
                            continue

                        interaction = int(section.split(" ")[0]) - 1

                        # Determinar la(s) clave(s) con las que trabajaremos
                        keys_to_update = []

                        if split_by_atom:
                            # Extraer los átomos del residuo: parte antes del primer '-'
                            # Formato: "1 | prot_atom1,prot_atom2-lig_atom |"
                            parts = section.split(GROUP_DELIM)
                            if len(parts) >= 2:
                                atoms_str = parts[1]  # "prot_atom1,prot_atom2-lig_atom"
                                prot_part = atoms_str.split("-")[0]  # "prot_atom1,prot_atom2"
                                atom_list = [a.strip() for a in prot_part.split(",") if a.strip()]
                                keys_to_update = [f"{residue_base} {atom}" for atom in atom_list]
                        else:
                            keys_to_update = [residue_base]

                        # Procesar cada clave (residuo o residuo+átomo)
                        for key in keys_to_update:
                            # Inicializar si no existe
                            if key not in data:
                                data[key] = nan_template.copy()

                            current_value = data[key][interaction]

                            # --- Lógica de acumulación (idéntica a la original) ---
                            if mode in ("min", "max") and activity > 0.0000:
                                if np.isnan(current_value) or op(activity, current_value):
                                    data[key][interaction] = activity
                            elif mode == "mean" and activity > 0.0000:
                                if not isinstance(current_value, list):
                                    data[key][interaction] = [1, activity]
                                else:
                                    count, total = current_value
                                    data[key][interaction] = [count + 1, total + activity]
                            elif mode in ("count", "percent"):
                                if np.isnan(current_value):
                                    data[key][interaction] = 1
                                else:
                                    data[key][interaction] += 1

            # Postprocesado para 'mean' y 'percent'
            if mode == "mean":
                for residue, interactions in data.items():
                    for i, value in enumerate(interactions):
                        if isinstance(value, list):
                            count, total = value
                            data[residue][i] = round(total / count, 2)
            elif mode == "percent":
                num_lines = len(matrix) - 1
                for residue, interactions in data.items():
                    for i, value in enumerate(interactions):
                        data[residue][i] = (value / num_lines * 100) if not np.isnan(value) else np.nan

            return data

        def initialize_data(matrix: list[list[str]]) -> tuple[list[float], list[list[str]]]:
            """
            Devuelve una plantilla de NaN (una por tipo de interacción) y la matriz transpuesta.

            Args:
                matrix (list[list[str]]): Matriz validada (filas = ligandos, columnas = residuos).

            Returns:
                tuple: (nan_template, transposed_matrix)
                    - nan_template: lista de np.nan con longitud igual al número de tipos de interacción.
                    - transposed_matrix: matriz transpuesta (filas = residuos, columnas = ligandos).
            """
            nan_template = [np.nan] * len(self.interaction_labels)
            transposed = self.transpose_matrix(interaction_data=matrix)
            return nan_template, transposed

        def plot_heatmap(self, data, title, x_label, y_label, mode, min_v, max_v, save, case, remove_empty):
            """
            Creates and optionally saves the heatmap visualization.

            Args:
                data (dict): The processed interaction data to visualize.
                title (str): The title for the heatmap.
                x_label (str): Label for the x-axis.
                y_label (str): Label for the y-axis.
                mode (str): The mode used to process the data.
                min_v (int, optional): Minimum value for the heatmap color scale.
                max_v (int, optional): Maximum value for the heatmap color scale.
                save (bool): Whether to save the heatmap to a file.

            Returns:
                None: Displays the heatmap or saves it to a file.
            """

            def residue_sort_key(label: str):
                """
                Clave de ordenación para los residuos (y opcionalmente átomos).
                Formato esperado: 'XXX num' o 'XXX num AA' (AA = átomo).
                Orden primario: número (segundo elemento tras split).
                Orden secundario: átomo (tercer elemento, alfabético; si no existe, cadena vacía).
                """
                parts = label.split()
                num = int(parts[1]) if len(parts) > 1 else 0
                atom = parts[2] if len(parts) > 2 else ""
                return (num, atom)

            df = pd.DataFrame(
                data,
                index=self.interaction_labels
                if case is None
                else [elemento.upper() for elemento in self.interaction_labels]
                if case == "upper"
                else [elemento.lower() for elemento in self.interaction_labels],
            )

            sorted_columns = sorted(df.columns, key=residue_sort_key)
            df = df.reindex(columns=sorted_columns)

            # Filtrar filas completamente vacías si se solicita
            if remove_empty:
                df = df.dropna(how="all")
                if df.empty:
                    print("Advertencia: No hay datos para mostrar después de eliminar las filas vacías.")
                    return

            max_cols = self.heat_max_cols  # Maximum number of columns per heatmap
            num_cols = len(df.columns)

            vmin = min_v if min_v else df.min().min()
            vmax = max_v if max_v else df.max().max()

            num_heatmaps = (num_cols + max_cols - 1) // max_cols  # Round up
            cols_per_heatmap = (num_cols + num_heatmaps - 1) // num_heatmaps  # Distribute evenly

            rang = float(vmax - vmin)

            if mode != "count":
                rang *= 10

            trunca = int(rang)
            ticks = min(6, trunca + 1)

            # if title ends with '(/)', it will not be displayed in the heatmap
            display = True
            if title.endswith("(/)"):
                title = title[:-3]
                display = False

            for i in range(num_heatmaps):
                start_col = i * cols_per_heatmap
                end_col = min((i + 1) * cols_per_heatmap, num_cols)
                df_subset = df.iloc[:, start_col:end_col]

                fig, ax = plt.subplots(figsize=(14, 9))
                sns.heatmap(
                    df_subset,
                    annot=True,
                    linewidths=0.5,
                    linecolor="lightgrey",
                    cmap=self.heat_colors,
                    fmt=".0f" if mode == "count" else ".1f",
                    vmin=vmin,
                    vmax=vmax,
                    cbar_kws={
                        "ticks": np.linspace(vmin, vmax, num=ticks),
                        "format": "%.0f" if mode == "count" else "%.2f",
                    },
                )

                fig.canvas.draw()
                # Add Subsite colors
                if residue_color_map:
                    for label in ax.get_xticklabels():
                        res = label.get_text().strip()

                        # Primero buscar coincidencia exacta
                        back_color = residue_color_map.get(res, "white")

                        # Si no hay coincidencia exacta, buscar versión con "+"
                        if back_color == "white":
                            # Verificar si el residuo tiene formato "RES NUM" (con dos partes)
                            parts = res.split()
                            if len(parts) >= 3:
                                # Probar con el formato "RES NUM +"
                                if parts[2] not in ["N", "CA", "C", "O"]:
                                    res_with_plus = f"{parts[0]} {parts[1]} +"
                                    back_color = residue_color_map.get(res_with_plus, "white")
                                if back_color == "white":
                                    res_with_plus = f"{parts[0]} {parts[1]}"
                                    back_color = residue_color_map.get(res_with_plus, "white")

                        # Si aún es blanco y el residuo es de cadena principal (tiene 3 partes: RES NUM ATOMO)
                        # En ese caso se queda blanco (como pide el requisito)

                        if back_color == "white":
                            label.set_color("black")
                        else:
                            label.set_color(self._get_contrasting_text_color(back_color))

                        label.set_bbox(dict(facecolor=back_color, edgecolor="none", boxstyle="round,pad=0.2"))
                        label.set_fontfamily("DejaVu Sans Mono")
                    if subsite_color_map:
                        legend_handles = [Patch(facecolor=subsite_color_map[s], label=s) for s in subsite_color_map]
                        fig.legend(
                            handles=legend_handles,
                            loc="upper center",
                            bbox_to_anchor=(0.5, 0.05),
                            ncol=min(5, len(legend_handles)),
                            frameon=False,
                        )

                # if title ends with '(/)', it will not be displayed in the heatmap
                if display:
                    if num_heatmaps == 1:
                        plt.title(f"{title}")
                    else:
                        plt.title(f"{title} (Columns {start_col + 1}-{end_col})")

                ax.set_xlabel(x_label)
                ax.set_ylabel(y_label)
                fig.subplots_adjust(left=0.145, bottom=0.155, right=1, top=0.935)

                if not save:
                    plt.show()
                else:
                    if num_heatmaps == 1:
                        filename = os.path.join(self.saving_directory, f"{title}.png")
                    else:
                        filename = os.path.join(self.saving_directory, f"{title}_part_{i + 1}.png")
                    plt.savefig(filename)
                    plt.close()

        self._verify_case(case=case)

        data = copy.deepcopy(interaction_data)
        matrix = data.matrix

        self.set_config(interaction_data=interaction_data)

        # Validate and prepare the matrix
        matrix = validate_and_prepare_matrix(matrix=matrix, mode=mode)

        # Ensure the mode is valid
        if mode not in HEATMAP_MODES:
            raise InvalidModeException(mode=mode, expected_values=HEATMAP_MODES)

        # Process the matrix based on the mode
        nan_template, matrix = initialize_data(matrix=matrix)
        data = process_matrix(matrix=matrix, mode=mode, nan_template=nan_template, split_by_atom=split_by_atom)

        # Generate and display/save the heatmap
        plot_heatmap(
            self=self,
            data=data,
            title=title,
            x_label=x_label,
            y_label=y_label,
            mode=mode,
            min_v=min_v,
            max_v=max_v,
            save=save,
            case=case,
            remove_empty=remove_empty,
        )

    def bar_chart(
        self,
        interaction_data: InteractionData,
        plot_name: str,
        axis: str = "rows",
        label_x: str = None,
        label_y: str = "Number of intermolecular interactions",
        title: str = "Protein-drug interactions",
        stacked: bool = False,
        save: bool = False,
        colors: list[str] = None,
        type_count: bool = False,
        subpocket_path: str = None,
        subpocket_colors: list[str] = None,
        export_xlsx: bool = False,
        case: str = None,
    ) -> None:
        """
        Generates a bar chart based on interaction data.

        This method extracts relevant data from an interaction matrix and visualizes it
        as a bar chart. If the dataset is too large, it automatically splits the data into
        multiple plots for better readability.

        The method supports:
        - **Standard bar charts**: A single bar per residue or PDB complex.
        - **Stacked bar charts** (`stacked=True`): Bars grouped by interaction types,
            showing their relative contribution.
        - **Interactive annotations**: When save=False, hovering over a bar displays
            context-specific details (e.g., absolute counts, percentages, or totals), depending
            on the plot type.
        - **Automatic data splitting**: Large datasets are split into multiple charts.

        Args:
            interaction_data (InteractionData): The object containing the interaction matrix.
            plot_name (str): Name of the plot (used for saving).
            axis (str): Defines whether to plot rows ('rows') or columns ('columns').
            label_x (str, optional): Label for the x-axis. Defaults to "Interacting protein residues".
            label_y (str, optional): Label for the y-axis. Defaults to "Number of intermolecular interactions".
            title (str, optional): Title of the chart. Defaults to "Protein-drug interactions".
            stacked (bool, optional): If True, creates a stacked bar chart. Defaults to False.
            save (bool, optional): If True, saves the chart as a PNG file. If False,
              the function displays an interactive plot. Interactive mode requires a
              GUI backend (e.g., by running %matplotlib qt or %matplotlib tk). Defaults to False.
            colors (list[str], optional): List of colors for interaction types. Defaults to None.
            type_count (bool, optional): Controls how interaction values are counted.
              If `True`, counts **all occurences** of each interaction type, even if repeated .
              If `False`, each interaction type is counted only **once per residue/complex**. Defaults to `False`.
            subpocket_path (str, optional): Path to a CSV defining subpockets and their residues.
               Used to color residue labels by subsite. If the file defines fewer than two
               subsites, coloring is disabled.
            subpocket_colors (list[str] or str, optional): Custom colors for subpockets
               (e.g., ['#8C2CE2','#DFBA52'], colors are used in order and repeated if needed),
               or a colormap name (e.g., "Dark2", "tab10". it is interpreted as a Matplotlib colormap name).
            export_xlsx (bool, optional): If True, exports the bar-chart data to an Excel file
               using the same base name as `plot_name` (with `.xlsx` extension). Defaults to False.
            case (str, optional): Case style for the plot leyend. Can be 'upper', 'lower', or None. Defaults to None.

        Returns:
            None: The function either displays or saves the plot.

        Raises:
            ValueError: If `axis` is not 'rows' or 'columns'.
        """

        self._verify_case(case=case)

        residue_color_map = None
        subsite_color_map = None

        if subpocket_path is not None:
            residue_color_map, subsite_color_map = self._load_subsite_color_maps(subpocket_path, subpocket_colors)
        matrix = interaction_data.matrix
        self.set_config(interaction_data=interaction_data)
        # Initialize and get data
        colors, data, indices, transposed_data = self._plot_init(colors, matrix, axis, type_count)
        if not stacked:
            y_values = [sum(row) for row in data]
            transposed_data = None
        else:
            y_values = None
        max_elements_plot = self.plot_max_cols
        num_x_elements = len(indices)

        # -------------------------------
        # Divide data into chunks if needed
        # -------------------------------
        plots = []
        if num_x_elements > max_elements_plot:
            for start in range(0, num_x_elements, max_elements_plot):
                end = start + max_elements_plot

                sub_indices = indices[start:end]

                if stacked:
                    sub_transposed = [group[start:end] for group in transposed_data]
                    plots.append((sub_indices, None, sub_transposed))
                else:
                    sub_y = y_values[start:end]
                    plots.append((sub_indices, sub_y, None))

        else:
            if stacked:
                plots = [(indices, None, transposed_data)]
            else:
                plots = [(indices, y_values, None)]

        for plot_i, (sub_indices, sub_data, sub_transposed) in enumerate(plots):
            # Original plotting logic if data fits in one plot
            fig, ax = plt.subplots(num=plot_name + f"_{plot_i}", figsize=(16, 6))
            if stacked:
                right_margin = 0.75
            else:
                right_margin = 0.95
            fig.subplots_adjust(left=0.10, right=right_margin, bottom=0.25, top=0.95)

            if stacked:
                bars = []
                bottoms = [0] * len(sub_indices)
                for index, group in enumerate(sub_transposed):
                    label = (
                        self.interaction_labels[index]
                        if case is None
                        else self.interaction_labels[index].upper()
                        if case == "upper"
                        else self.interaction_labels[index].lower()
                    )
                    bars.append(ax.bar(sub_indices, group, bottom=bottoms, label=label, color=colors[index]))
                    bottoms = [i + j for i, j in zip(bottoms, group)]
                max_y = max([sum(col) for col in zip(*sub_transposed)]) if sub_transposed and sub_transposed[0] else 0
                ax.legend(loc="center left", bbox_to_anchor=(1, 0.5), ncol=1)
                # Add interactive tooltips only when showing figure (not when saving)
                if not save:
                    cursor = mplcursors.cursor(bars, hover=True)

                    @cursor.connect("add")
                    def on_add(sel):
                        index = sel.index
                        # Absolute counts
                        abs_counts = [sub_transposed[i][index] for i in range(len(sub_transposed))]
                        total = sum(abs_counts)
                        # Compute percentages using ONLY the current chunk (sub_transposed), not the whole data
                        percentages = [
                            (abs_counts[i] / total * 100) if total > 0 else 0 for i in range(len(abs_counts))
                        ]
                        annotation_lines = [f"TOTAL interactions: {total}"]
                        annotation_lines += [
                            f"{self.interaction_labels[i]}: {abs_counts[i]} ({percentages[i]:.1f} %)"
                            for i in range(len(abs_counts))
                        ]
                        sel.annotation.set_text("\n".join(annotation_lines))
                        sel.annotation.get_bbox_patch().set(fc="white", alpha=0.9)

            else:
                ax.bar(sub_indices, sub_data, color=colors[0] if colors else None)
                max_y = max(sub_data) if sub_data else 0
                if not save:
                    cursor = mplcursors.cursor(ax.containers[0], hover=True)

                    @cursor.connect("add")
                    def on_add(sel):
                        idx = sel.index
                        x_lab = sub_indices[idx]
                        y_val = sub_data[idx]
                        msg = (
                            f"{x_lab}: {y_val} total interactions"
                            if type_count
                            else f"{x_lab}: {y_val} distinct-type interactions"
                        )
                        sel.annotation.set_text(msg)
                        sel.annotation.get_bbox_patch().set(fc="white", alpha=0.9)

            ax.set_ylim(0, (max_y * 1.1) if max_y > 0 else 1)
            ax.tick_params(axis="x", labelrotation=90)
            # Add subsites colors
            if residue_color_map is not None and any(idx.split("-")[0] in residue_color_map for idx in sub_indices):
                for label in ax.get_xticklabels():
                    residue_name = label.get_text()
                    cleaned = residue_name.split("-")[0]  # Remove chain ("-A") if present
                    back_color = residue_color_map.get(cleaned, "white")
                    if back_color == "white":
                        label.set_color("black")
                    else:
                        label.set_color(self._get_contrasting_text_color(back_color))
                    label.set_bbox(dict(facecolor=back_color, edgecolor="none", boxstyle="round,pad=0.2"))
                    label.set_fontfamily("DejaVu Sans Mono")
                if subsite_color_map:
                    legend_handles = [Patch(facecolor=subsite_color_map[s], label=s) for s in subsite_color_map]
                    fig.legend(
                        handles=legend_handles,
                        loc="upper center",
                        bbox_to_anchor=(0.5, 0.1),
                        ncol=min(5, len(legend_handles)),
                        frameon=False,
                    )
            ax.set_ylabel(label_y)
            ax.set_xlabel(label_x or "Interacting protein residues")
            ax.set_title(title)
            ax.yaxis.set_major_locator(MaxNLocator(integer=True))
            if len(plots) == 1:
                out_name = plot_name
            else:
                out_name = f"{plot_name}_part{plot_i + 1}"
            self._plot_end(save, plt, fig, out_name)
        # Export data to Excel if requested
        if export_xlsx:
            xlsx_name = f"{plot_name}.xlsx"
            xlsx_path = os.path.join(self.saving_directory, xlsx_name)

            if not stacked:
                self._export_bar_data_to_excel(
                    filename=xlsx_path,
                    indices=indices,
                    data=[[v] for v in y_values],
                    interaction_labels=["Total_raw_interactions"]
                    if type_count
                    else ["Total_distinct_interaction_types"],
                )

            # CASE 2 — stacked=True (regardless of type_count)
            elif stacked:
                self._export_bar_data_to_excel(
                    filename=xlsx_path, indices=indices, data=data, interaction_labels=self.interaction_labels
                )

    def pie_chart(
        self,
        interaction_data: InteractionData,
        plot_name: str,
        axis: str,
        title: str = "Interaction Types",
        save: bool = False,
        colors: list[str] = None,
        type_count: bool = False,
        case: str = None,
        export_xlsx: bool = False,
    ) -> None:
        """
        Generates a pie chart based on interaction data.

        This method processes an interaction matrix, calculates the proportion of different interaction
        types, and visualizes them as a pie chart. Interaction types with zero occurrences are automatically
        filtered out. The chart can be displayed or saved as a PNG file.

        Args:
            interaction_data (InteractionData): The object containing the interaction matrix.
            plot_name (str): The name of the plot (used for saving).
            axis (str): Defines whether to analyze rows ('rows') or columns ('columns').
            title (str): Title plot
            save (bool, optional): If True, saves the pie chart instead of displaying it.
                                   If False, an interactive chart is geenrated. Defaults to False.
            colors (list[str], optional): List of colors for interaction types. Defaults to None.
            type_count (bool, optional): If False, counts distinct interaction types
                                        If True, sums all interaction values. Defaults to False.
            case (str, optional): Case style for the plot leyend. Can be 'upper', 'lower', or None. Defaults to None.
            export_xlsx (`bool`, optional): If True, exports the bar-chart data to an Excel file using the same name
                                            defined in `plot_name` (with `.xlsx` extension). Defaults to `False`.

        Returns:
            None: The function either displays or saves the plot.
        """

        def _filter_non_zero_interactions(total_interactions: list[int], colors: list[str]) -> list[tuple]:
            """
            Filters out interaction types with zero total count.

            Args:
                total_interactions (list[int]): Total interactions per type.
                colors (list[str]): List of colors for each interaction type.

            Returns:
                list[tuple]: Filtered interaction data as (label, total, color).
            """
            return [
                (label, total, colors[i])
                for i, (label, total) in enumerate(zip(self.interaction_labels, total_interactions))
                if total > 0
            ]

        def _plot_pie_chart(
            labels: list[str], sizes: list[int], colors: list[str], total_interactions: list[int], case: str, title: str
        ) -> None:
            """
            Creates and formats a pie chart.

            Args:
                labels (list[str]): Interaction labels for the pie chart.
                sizes (list[int]): Interaction sizes (counts or percentages).
                colors (list[str]): Colors for each interaction type.
                total_interactions (list[int]): Total counts of interactions per type.

            Returns:
                None
            """
            fig, ax_pie = plt.subplots(figsize=(10, 6))

            # Plot the pie chart
            wedges, texts = ax_pie.pie(sizes, labels=None, colors=colors, startangle=140)
            ax_pie.set_title(title)
            for w in wedges:
                w.set_picker(5)

            # Calculate percentages and create legend
            total = sum(total_interactions)
            legend_labels = [
                f"{label if case is None else label.upper() if case == 'upper' else label.lower()} "
                f"({round(count / total * 100, 2)}%)"
                for label, count in zip(self.interaction_labels, total_interactions)
                if count != 0
            ]
            ax_pie.legend(legend_labels, title="Interaction Types", loc="center left", bbox_to_anchor=(1, 0, 0.5, 1))

            wedges = list(wedges)

            return fig, ax_pie, wedges

        def _prepare_pie_chart_data(non_zero_interactions: list[tuple]) -> tuple:
            """
            Prepares data for plotting the pie chart.

            Args:
                non_zero_interactions (list[tuple]): Interaction data with non-zero totals.

            Returns:
                tuple: Labels, sizes, and colors for the pie chart.
            """
            if non_zero_interactions:
                return zip(*non_zero_interactions)  # Separates into labels, sizes, and colors
            return [], [], []

        self._verify_case(case=case)
        self.set_config(interaction_data=interaction_data)
        matrix = interaction_data.matrix

        # Initialize plotting parameters and data
        colors, data, indices, transposed_data = self._plot_init(colors, matrix, axis, type_count)

        # Calculate total interaction and prepade data for plotting
        total_interactions = [sum(transposed_data[i]) for i in range(len(transposed_data))]
        non_zero_interactions = _filter_non_zero_interactions(total_interactions, colors)
        labels_pie, sizes, pie_colors = _prepare_pie_chart_data(non_zero_interactions)
        labels_pie = list(labels_pie)
        sizes = list(sizes)
        pie_colors = list(pie_colors)

        if not sizes:
            print("No interaction types with non-zero counts to plot.")
            return
        # Plot the pie chart
        fig, ax_pie, wedges = _plot_pie_chart(labels_pie, sizes, pie_colors, total_interactions, case, title)

        if (not save) and wedges:
            # Create one persistent annotation (hidden by default)
            annot = ax_pie.annotate(
                "",
                xy=(0, 0),
                xytext=(10, 10),
                textcoords="offset points",
                bbox=dict(boxstyle="round", fc="white", alpha=0.9),
                arrowprops=dict(arrowstyle="->"),
            )
            annot.set_visible(False)

            def _on_move(event):
                if event.inaxes != ax_pie:
                    if annot.get_visible():
                        annot.set_visible(False)
                        fig.canvas.draw_idle()
                    return

                for i, w in enumerate(wedges):
                    contains, _ = w.contains(event)
                    if contains:
                        label = labels_pie[i]
                        count = sizes[i]
                        total = sum(sizes)
                        pct = (count / total * 100) if total else 0

                        annot.xy = (event.xdata, event.ydata)
                        annot.set_text(f"{label}: {count} ({pct:.1f}%)")
                        annot.set_visible(True)
                        fig.canvas.draw_idle()
                        return

                if annot.get_visible():
                    annot.set_visible(False)
                    fig.canvas.draw_idle()

            cid = fig.canvas.mpl_connect("motion_notify_event", _on_move)

            # Keep references to avoid garbage collection
            if not hasattr(self, "_pie_hover_state"):
                self._pie_hover_state = {}
            self._pie_hover_state[plot_name] = {"cid": cid, "annot": annot, "fig": fig}
        # Show or save the plot
        self._plot_end(save, plt, fig, plot_name)

        if export_xlsx:
            xlsx_path = os.path.join(self.saving_directory, f"{plot_name}.xlsx")
            self._export_pie_data_to_excel(xlsx_path, labels_pie, sizes)
