"""Filtering, sorting and shaping operations on interaction matrices.

Extracted verbatim from the original monolithic ``analyze_interactions.py``
(Fase 2, paso 4 del plan de modularización). No logic changes were made
during this extraction — only the physical location of the code changed.

Includes ``_get_interactions`` and ``_stack_reactives`` here rather than in
``_validation.py``: per Anexo A these two were flagged as "compartido,
revisar uso real al extraer" — checking actual usage showed
``_get_interactions`` is only ever called from ``_stack_reactives``, and
Anexo A itself places ``_stack_reactives`` in ``filter_mixin.py``. So both
land here together, next to their only caller chain.
"""

import copy
import csv
import os
import re

from .constants import (
    DIFF_DELIM,
    EMPTY_DASH_CELL,
    GROUP_DELIM,
    SAME_DELIM,
    is_not_empty_or_dash,
)
from .exceptions import InvalidAxisException, InvalidModeException, MissedActivityException
from .models import InteractionData


class FilterMixin:
    """Filtering, sorting, reshaping and interaction-counting operations.

    Depends on ``ValidationMixin`` (``self._check_variable_types``,
    ``self._verify_dimensions``, ``self._get_residues_axis``) and on
    ``self.save_interaction_data`` (still defined directly on
    ``AnalyzeInteractions`` at this point in the modularization; will move
    to ``export_mixin.py`` in a later step). Also reads ``self.interaction_labels``.
    """

    def _get_interactions(self, cell: str) -> list[int]:
        """
        Extracts and counts interactions from a cell string.

        Args:
            cell (str): The cell string containing interaction data.

        Returns:
            list: A list of interaction counts for each interaction type.
        """
        interactions = [0] * len(self.interaction_labels)
        sections = cell.split(GROUP_DELIM)
        for index in range(1, len(sections), 2):
            interaction = int(sections[index - 1].replace(DIFF_DELIM, "").replace(" ", ""))
            interactions[interaction - 1] += len(sections[index].split(SAME_DELIM))
        return interactions

    def _stack_reactives(
        self, matrix: list[list[str]], axis: str, type_count: bool
    ) -> tuple[list[list[int]], list[str]]:
        """
        Computes interaction-type counts for each element (residue or PDB entry),
        either by rows or by columns depending on `axis`.

        - If `axis='rows'`, each matrix row corresponds to an element.
        - If `axis='columns'`, the matrix is transposed so that each original column
          becomes an element, and interactions are accumulated per column instead of per row.

        The function returns:
            - A list of lists, where each sublist contains the counts of each
              interaction type for one element.
            - A list of element labels (indices), aligned with the rows of the output matrix.

        Args:
            matrix (list of lists): The interaction matrix.
            axis (str): 'rows' to process residues; 'columns' to process complexes.
            type_count (bool):
                - If True, counts *all occurrences* of each interaction type.
                - If False, counts *presence/absence* of each interaction type per cell.

        Returns:
            tuple[list[list[int]], list[str]]:
                (interaction_counts, element_labels)
        """
        self._verify_dimensions(matrix=matrix)

        # If we want columns, transpose so that rows become the elements to count
        if axis == "columns":
            matrix = self.transpose_matrix(matrix)

        # Number of elements (skip header row)
        num_elements = len(matrix) - 1
        num_types = len(self.interaction_labels)

        # Initialize result structure
        reactives = [[0] * num_types for _ in range(num_elements)]
        indices = [matrix[row][0].split("_")[0].strip() for row in range(1, len(matrix))]

        # Count interactions
        for row in range(1, len(matrix)):
            for col in range(1, len(matrix[row])):
                cell = matrix[row][col]
                interactions = self._get_interactions(cell)

                for i in range(num_types):
                    if type_count:
                        reactives[row - 1][i] += interactions[i]
                    elif interactions[i] > 0:
                        reactives[row - 1][i] += 1

        return reactives, indices

    def filter_by_interaction(
        self, interaction_data: InteractionData, interactions: list[int], save: str = None
    ) -> InteractionData:
        """
        Filters an interaction matrix based on specified interaction types.

        Args:
            interaction_data (InteractionData): The object containing the interaction matrix.
            interactions (list[int]): List of valid interaction types (numbers 1 to 7) to retain in the matrix.
            save (str, optional): File path to save the filtered matrix. Defaults to None.

        Returns:
            InteractionData: The updated InteractionData object with the filtered matrix.

        Raises:
            ValueError: If the matrix dimensions are invalid, or if no matching interactions are found.
        """

        def validate_list(interactions: list[int], interaction_labels: list[str]) -> None:
            """
            Validates the interaction list to ensure it contains unique numbers between 1 and 7.

            Args:
                interactions (list): List of integers representing interaction types.

            Raises:
                ValueError: If any number is outside the range of 1 to 7 or if there are duplicates.
            """
            # Valid interactions are numbers from 1 to 7
            valid_numbers = set(range(1, len(interaction_labels) + 1))

            # Check if all numbers in the list are within the valid range
            for num in interactions:
                if num not in valid_numbers:
                    raise ValueError(
                        f"Invalid interaction: {num}. Must be a number between "
                        f"{min(valid_numbers)} and {max(valid_numbers)}."
                    )

            # Ensure the list contains no duplicate values
            if len(set(interactions)) != len(interactions):
                raise ValueError("The interaction list must not contain duplicates.")

        # Validate types of the matrix, interactions, and save parameters
        self._check_variable_types(
            variables=[interaction_data, interactions, save],
            expected_types=[InteractionData, list, (str, None.__class__)],
            variable_names=["interaction_data", "interactions", "save"],
        )

        data = copy.deepcopy(interaction_data)
        matrix = data.matrix

        # Validate that the matrix has appropriate dimensions
        self._verify_dimensions(matrix=matrix)

        # Validate that the interaction list contains valid values
        validate_list(interactions=interactions, interaction_labels=data.interactions)

        # Track whether any interactions were filtered
        changes = False

        # Iterate through each cell in the matrix (skipping the header row/column)
        for i in range(1, len(matrix)):
            for j in range(1, len(matrix[i])):
                cell = matrix[i][j]

                # If the cell is not empty ('-'), process it
                if is_not_empty_or_dash(cell):
                    sections = cell.split(DIFF_DELIM)
                    cell = ""

                    # Iterate through the sections in the cell to filter valid interactions
                    for section in sections:
                        # Check if the first number in the section is in the valid interaction list
                        if int(section.split(" ")[0]) in interactions:
                            # Add the section to the cell if it contains a valid interaction
                            if cell == "":
                                cell = section
                            else:
                                cell += DIFF_DELIM + section
                            changes = True

                    # Update the cell with the filtered sections or set it to '-' if empty
                    matrix[i][j] = cell if cell != "" else EMPTY_DASH_CELL

        # If no changes were made, raise an error indicating no matching interactions were found
        if not changes:
            raise ValueError("No matching interactions were found in the matrix.")

        data.matrix = matrix

        # Save the filtered matrix to a file if a save path is provided
        if save:
            self.save_interaction_data(interaction_data=data, filename=save)

        return data

    def filter_by_subunit(
        self, interaction_data: InteractionData, subunits: list[str], save: str = None
    ) -> InteractionData:
        """
        Filters an interaction matrix based on specified subunits.

        This method processes an interaction matrix and removes rows or interaction elements
        that do not match the provided subunit list. The filtering mechanism depends on whether
        the matrix is organized by residues or interactions.

        Args:
            interaction_data (InteractionData): The object containing the interaction matrix.
            subunits (list[str]): List of valid subunits used as filtering criteria.
            save (str, optional): File path to save the filtered matrix. Defaults to None.

        Returns:
            InteractionData: The updated InteractionData object with the filtered matrix.

        Raises:
            ValueError: If the matrix dimensions are invalid or if no matching subunits are found.
        """

        def get_subunits_location(matrix: list[list[str]]) -> str:
            """
            Determines the location of subunits in the matrix (residues or interactions).

            Args:
                matrix (list): The matrix being analyzed.

            Returns:
                str: 'residues' if the first column indicates residue data, otherwise 'interactions'.
            """
            return "residues" if len(matrix[1][0].split("-")) == 2 else "interactions"

        # Check types of the matrix, subunits, and save parameters
        self._check_variable_types(
            variables=[interaction_data, subunits, save],
            expected_types=[InteractionData, list, (str, None.__class__)],
            variable_names=["interaction_data", "subunits", "save"],
        )

        data = copy.deepcopy(interaction_data)
        matrix = data.matrix

        # Validate the dimensions of the matrix
        self._verify_dimensions(matrix=matrix)

        # Determine the axis of residues in the matrix
        axis = self._get_residues_axis(matrix=matrix)

        # Transpose the matrix if the axis is columns
        if axis == "columns":
            matrix = self.transpose_matrix(matrix)

        # Determine whether the matrix contains residues or interactions
        subunitsLocation = get_subunits_location(matrix=matrix)

        # Initialize change tracking variables
        changes = 0

        # Filter based on residue locations
        if subunitsLocation == "residues":
            for index in range(1, len(matrix)):
                sections = matrix[index - changes][0].split("-")
                # Remove rows without valid sections or subunits
                if len(sections) != 2 or sections[1] not in subunits:
                    matrix.pop(index - changes)
                    changes += 1
        else:
            # Iterate through each cell in the matrix to filter interactions
            for i in range(1, len(matrix)):
                for j in range(1, len(matrix[i])):
                    cell = matrix[i][j]

                    # Process non-empty cells
                    if is_not_empty_or_dash(cell=cell):
                        sections = cell.split(DIFF_DELIM)
                        cell = ""

                        # Iterate through each section in the cell
                        for section in sections:
                            separators = section.split(GROUP_DELIM)[:-1]

                            # Filter out unwanted interactions
                            for index in range(1, len(separators)):
                                if index % 2 != 0:
                                    interactions = separators[index].split(SAME_DELIM)
                                    subchanges = 0

                                    # Remove interactions not in the valid subunits
                                    for interaction in range(len(interactions)):
                                        raw = interactions[interaction - subchanges]
                                        if "(" in raw and ")" in raw:
                                            subunit = raw.partition("(")[2].rpartition(")")[0]
                                        else:
                                            subunit = None
                                        if subunit not in subunits:
                                            changes += 1
                                            interactions.pop(interaction - subchanges)
                                            subchanges += 1

                                    # Rebuild the cell if there are valid interactions
                                    if interactions:
                                        cell += separators[index - 1] + GROUP_DELIM
                                        cell += SAME_DELIM.join(interactions) + GROUP_DELIM + DIFF_DELIM

                        # Update the cell with filtered interactions
                        cell = cell[:-2]  # Remove trailing comma and space
                        matrix[i][j] = cell if cell else EMPTY_DASH_CELL

        # Transpose the matrix back if it was originally in columns
        if axis == "columns":
            matrix = self.transpose_matrix(matrix)

        data.matrix = matrix
        data.subunits_set = subunits

        # Save the filtered matrix to a file if a save path is provided
        if save:
            self.save_interaction_data(interaction_data=data, filename=save)

        return data

    def filter_by_residue(
        self,
        interaction_data: InteractionData,
        chain: str = None,
        subpocket_path: str = None,
        subpockets: list[str] = None,
        save: str = None,
    ) -> InteractionData:
        """
        Filters an interaction matrix based on a specified residue or subpockets.

        This method processes an interaction matrix and removes rows or interaction elements
        that do not match the provided chain or residues extracted from the given subpockets.

        The filtering can be performed in two ways:
        1. By specifying a `chain` ("<main>" or "<side>"), which filters interactions
        based on the presence of main or side chain atoms.
        2. By providing a `subpocket_path` and a list of `subpockets`, which extracts
        residues from a predefined subpocket file and filters interactions accordingly.

        Args:
            interaction_data (InteractionData): The object containing the interaction matrix.
            chain (str, optional): Specifies whether to retain "<main>" or "<side>" interactions. Defaults to None.
            subpocket_path (str, optional): Path to the file containing subpocket residue definitions. Defaults to None.
            subpockets (list[str], optional): List of subpockets to use for residue-based filtering. Defaults to None.
            save (str, optional): File path to save the filtered matrix. Defaults to None.

        Returns:
            InteractionData: The updated InteractionData object with the filtered matrix.

        Raises:
            ValueError: If the matrix dimensions are invalid.
            InvalidModeException: If an invalid chain mode is specified.
            Exception: If no protein atoms are available in the interaction data when filtering by chain.
        """
        MAIN_ATOMS = {"C", "CA", "N", "O"}  # átomos de cadena principal

        def extract_subpockets_from_file(subpocket_file_path: str, subpocket_list: list[str]) -> list[dict]:
            """
            Extrae residuos del fichero de subpockets, devolviendo una lista de diccionarios
            con la información de residuo, átomo y cadena.

            Formato esperado del fichero (una fila por subpocket): el nombre del
            subpocket en la primera columna y los residuos que lo componen en las
            columnas siguientes, separados por comas (p. ej. "S1,ARG32,MET54").
            Los espacios en blanco alrededor de cada valor se ignoran.
            """
            residues = []
            with open(subpocket_file_path, mode="r", encoding="utf-8") as csv_file:
                reader = csv.reader(csv_file)
                for row in reader:
                    if not row:
                        continue
                    if row[0].strip() in subpocket_list:
                        items = [item.strip() for item in row[1:]]
                        for item in items:
                            if not item:
                                continue
                            # 1. Extraer cadena (<main>, <side>) si existe
                            chain = None
                            if "<" in item:
                                base, chain_suffix = item.split("<", 1)
                                chain = "<" + chain_suffix
                            else:
                                base = item
                            # 2. Separar residuo y posible átomo
                            parts = base.split()
                            if len(parts) == 2:
                                residue_full, atom = parts[0], parts[1]
                            elif len(parts) == 1:
                                residue_full, atom = parts[0], None
                            else:
                                continue  # formato no esperado
                            # 3. Extraer nombre y número de residuo (ej: THR25)
                            match = re.match(r"([A-Za-z]+)(\d+)", residue_full)
                            if not match:
                                continue
                            residue_id = match.group(1) + match.group(2)  # "THR25"
                            residues.append({"residue_id": residue_id, "atom": atom, "chain": chain})
            return residues

        def filter_cell(cell: str, atom_set: set, chain: str) -> str:
            """
            Filtra una celda conservando solo los pares de interacción cuyos átomos
            de la proteína cumplan los criterios de átomo o cadena.
            """
            if not is_not_empty_or_dash(cell):
                return cell
            new_parts = []
            interactions = cell.split(DIFF_DELIM)
            for interaction in interactions:
                if not interaction.strip():
                    continue
                parts = interaction.split(GROUP_DELIM, 1)
                if len(parts) < 2:
                    continue
                interaction_number = parts[0]
                pairs_str = parts[1]
                valid_pairs = []
                for pair in pairs_str.split(SAME_DELIM):
                    prot_part = pair.split("-")[0]
                    prot_atoms = [a.strip() for a in prot_part.split(",")]
                    keep = False
                    for atom in prot_atoms:
                        if atom_set is not None and atom in atom_set:
                            keep = True
                            break
                        if chain is not None:
                            is_main = atom in MAIN_ATOMS
                            if (chain == "<main>" and is_main) or (chain == "<side>" and not is_main):
                                keep = True
                                break
                    if keep:
                        valid_pairs.append(pair)
                if valid_pairs:
                    new_parts.append(f"{interaction_number}{GROUP_DELIM}{SAME_DELIM.join(valid_pairs)}")
            if new_parts:
                return DIFF_DELIM.join(new_parts)
            return "-"

        def filter_matrix(matrix: list[list[str]], residues: list[dict], global_chain: str) -> list[list[str]]:
            """
            Filtra la matriz conservando únicamente las filas y los elementos de interacción
            que coinciden con los residuos (y sus átomos/cadenas) extraídos de los subpockets.
            """
            # Construir las restricciones para cada residuo
            constraints = {}  # residue_id -> {'keep_all': bool, 'atoms': set or None, 'chain': str or None}
            for entry in residues:
                rid = entry["residue_id"]
                if rid not in constraints:
                    constraints[rid] = {"keep_all": False, "atoms": set(), "chain": None}

                # Si la entrada no especifica átomo ni cadena → mantener todas las interacciones del residuo
                if entry["atom"] is None and entry["chain"] is None:
                    constraints[rid]["keep_all"] = True
                    constraints[rid]["atoms"] = None
                    constraints[rid]["chain"] = None
                    continue  # las siguientes entradas para este residuo son irrelevantes

                # Si ya decidimos mantener todas, ignoramos cualquier otra restricción para este residuo
                if constraints[rid]["keep_all"]:
                    continue

                # Registrar átomo concreto, si existe
                if entry["atom"] is not None:
                    constraints[rid]["atoms"].add(entry["atom"])

                # Registrar cadena, si existe (asumimos consistencia; si no, prevalece la primera)
                if entry["chain"] is not None and constraints[rid]["chain"] is None:
                    constraints[rid]["chain"] = entry["chain"]

            filtered_matrix = [matrix[0]]  # conservar cabecera
            for row in matrix[1:]:
                residue_full = row[0].replace(" ", "").split("-")[0]
                # Obtener identificador base del residuo (ej. THR25)
                match = re.match(r"([A-Za-z]+\d+)", residue_full)
                res_id = match.group(1) if match else residue_full

                if res_id in constraints:
                    constr = constraints[res_id]
                    # Caso 1: mantener toda la fila
                    if constr["keep_all"]:
                        filtered_matrix.append(row)
                        continue

                    # Determinar criterios de filtrado
                    atom_set = constr["atoms"] if constr["atoms"] else None  # None si vacío
                    chain_filt = (
                        constr["chain"] if not atom_set else None
                    )  # átomos concretos tienen prioridad sobre cadena

                    new_row = [row[0]]
                    row_empty = True
                    for cell in row[1:]:
                        filtered_cell = filter_cell(cell, atom_set, chain_filt)
                        new_row.append(filtered_cell)
                        if filtered_cell != "-":
                            row_empty = False
                    if not row_empty:
                        filtered_matrix.append(new_row)

                elif global_chain:
                    # Aplicar filtro global de cadena a filas no incluidas en subpockets
                    new_row = [row[0]]
                    row_empty = True
                    for cell in row[1:]:
                        filtered_cell = filter_cell(cell, None, global_chain)
                        new_row.append(filtered_cell)
                        if filtered_cell != "-":
                            row_empty = False
                    if not row_empty:
                        filtered_matrix.append(new_row)
                # else: descartar la fila
            return filtered_matrix

        def validate_chain(chain: str) -> None:
            """Valida que la cadena sea '<main>' o '<side>'."""
            valid_chains = ["<main>", "<side>"]
            if chain not in valid_chains:
                raise InvalidModeException(mode=chain, expected_values=valid_chains)

        # Check types of the matrix, chain, and subpocket
        self._check_variable_types(
            variables=[interaction_data, chain, subpockets, subpocket_path, save],
            expected_types=[
                InteractionData,
                (str, None.__class__),
                (list, None.__class__),
                (str, None.__class__),
                (str, None.__class__),
            ],
            variable_names=["interaction_data", "chain", "subpockets", "subpocket_path", "save"],
        )

        filtered_data = copy.deepcopy(interaction_data)
        matrix = filtered_data.matrix
        residues_selection = []

        # Extraer residuos si se proporcionan subpockets
        if subpocket_path and subpockets:
            subpocket_path = os.path.join(self.input_directory, subpocket_path)
            residues = extract_subpockets_from_file(subpocket_path, subpockets)

            if chain:
                validate_chain(chain=chain)
                if not filtered_data.protein:
                    raise Exception("No protein atoms available in the interaction data.")
                # Filtrar las entradas según la cadena externa
                filtered_residues = []
                for entry in residues:
                    if entry["atom"] is not None:
                        # el átomo determina si es compatible con la cadena
                        atom_is_main = entry["atom"] in MAIN_ATOMS
                        if (chain == "<main>" and atom_is_main) or (chain == "<side>" and not atom_is_main):
                            filtered_residues.append(entry)
                    elif entry["chain"] is not None:
                        if entry["chain"] == chain:
                            filtered_residues.append(entry)
                    else:
                        # sin átomo ni cadena -> forzar la cadena externa
                        filtered_residues.append({"residue_id": entry["residue_id"], "atom": None, "chain": chain})
                residues_selection = filtered_residues
                chain = None  # ya se aplicó en la selección
            else:
                residues_selection = residues
        elif chain:
            validate_chain(chain=chain)
            if not filtered_data.protein:
                raise Exception("No protein atoms available in the interaction data.")

        # Validar dimensiones de la matriz
        self._verify_dimensions(matrix=matrix)

        axis = self._get_residues_axis(matrix=matrix)

        if axis == "columns":
            matrix = self.transpose_matrix(matrix)

        filtered_matrix = filter_matrix(matrix=matrix, residues=residues_selection, global_chain=chain)

        if axis == "columns":
            filtered_matrix = self.transpose_matrix(filtered_matrix)

        filtered_data.matrix = filtered_matrix

        if save:
            self.save_interaction_data(interaction_data=filtered_data, filename=save)

        return filtered_data

    def remove_empty_axis(self, interaction_data: InteractionData, save: str = None) -> InteractionData:
        """
        Removes empty rows and columns from the interaction matrix.

        This method iterates through the interaction matrix and removes rows and columns
        that contain only empty or placeholder values (e.g., dashes or empty strings).
        The cleaned matrix maintains the original structure and is optionally saved to a file.

        Args:
            interaction_data (InteractionData): The input interaction data from which empty
                rows and columns will be removed.
            save (str, optional): The filename to save the cleaned matrix. If None, the matrix will not be saved.

        Returns:
            InteractionData: The cleaned interaction data with empty rows and columns removed.

        Raises:
            TypeMismatchException: If the types of input variables do not match the expected types.
            ValueError: If the matrix dimensions are invalid or if it is too small after cleaning.
        """

        def _remove_empty_rows(matrix: list[list[str]]) -> list[list[str]]:
            """
            Helper function to remove empty rows from the matrix.

            Args:
                matrix (list[list[str]]): The matrix from which empty rows will be removed.

            Returns:
                list[list[str]]: The matrix with empty rows removed.
            """
            changes = 0
            for row in range(1, len(matrix)):
                if all(not is_not_empty_or_dash(cell=column) for column in matrix[row - changes][1:]):
                    matrix.pop(row - changes)
                    changes += 1
            return matrix

        self._check_variable_types(
            variables=[interaction_data, save],
            expected_types=[InteractionData, (str, None.__class__)],
            variable_names=["interaction_data", "save"],
        )

        data = copy.deepcopy(interaction_data)
        matrix = data.matrix

        self._verify_dimensions(matrix=matrix)

        matrix = _remove_empty_rows(matrix=matrix)

        # Transpose the matrix, remove empty columns (which are now rows)
        matrix = self.transpose_matrix(matrix)
        matrix = _remove_empty_rows(matrix=matrix)

        # Transpose back to restore original format
        matrix = self.transpose_matrix(matrix)

        data.matrix = matrix

        if save:
            self.save_interaction_data(interaction_data=data, filename=save)

        return data

    def sort_matrix(
        self,
        interaction_data: InteractionData,
        axis: str = "rows",
        thr_interactions: int = None,
        thr_activity: float = None,
        selected_items: int = None,
        _count: bool = False,
        residue_chain: bool = False,
        save: str = None,
    ) -> InteractionData:
        """
        Sorts and filters rows or columns in the interaction matrix based on interaction criteria.

        This method sorts and selects reactive rows or columns in the interaction matrix according
        to the specified criteria. It supports selection based on:
        - **Minimum interaction count** (`thr_interactions`): Keeps rows/columns with at least
        the specified number of interactions.
        - **Activity threshold** (`thr_activity`): Selects rows/columns where the activity value
        meets or exceeds the given threshold.
        - **Top-ranked selection** (`selected_items`): Retains only the top N rows/columns based
        on interaction count.
        - **Residue chain sorting** (`residue_chain=True`): Orders the matrix by residue index
        after filtering.

        The method can also return just the interaction count per row/column (`count=True`). If
        both `thr_interactions` and `selected_items` are provided, an error is raised.

        Args:
            interaction_data (InteractionData): The interaction data to be sorted.
            axis (str, optional): Specifies whether to sort rows ('rows') or columns ('columns'). Defaults to 'rows'.
            thr_interactions (int, optional): Minimum number of interactions required to retain a row/column.
            thr_activity (float, optional): Minimum activity value required to retain a row/column.
            selected_items (int, optional): Number of top rows/columns to keep based on interaction count.
            _count (bool, optional): If True, returns the count of interactions instead of modifying the matrix.
            residue_chain (bool, optional): If True, sorts the resulting matrix based on residue order in the chain.
            save (str, optional): File path to save the resulting matrix. Defaults to None.

        Returns:
            InteractionData: The sorted interaction matrix.

        Raises:
            ValueError: If multiple selection criteria (`thr_interactions`, `thr_activity`,
                `selected_items`) are used simultaneously.
            InvalidAxisException: If an invalid axis is provided.
            ValueError: If the matrix dimensions are insufficient for sorting.
        """

        def get_interactions(cell: str) -> int:
            """
            Counts the number of interactions in a cell formatted with interaction data.

            Args:
                cell (str): Cell containing interaction data, formatted with '|' separating values.

            Returns:
                int: Total number of interactions in the cell.
            """
            interactions = 0
            sections = cell.split(GROUP_DELIM)
            for index in range(1, len(sections), 2):
                interactions += len(sections[index].split(SAME_DELIM))
            return interactions

        def sort_by_residue(matrix: list[list[str]]) -> list[list[str]]:
            """
            Sorts the matrix based on residue indices in the first column.

            Args:
                matrix (list[list[str]]): The matrix to be sorted.

            Returns:
                list[list[str]]: The matrix sorted by residue indices.
            """
            # Validate matrix dimensions
            self._verify_dimensions(matrix=matrix)

            # Determine whether to sort by rows or columns
            axis = self._get_residues_axis(matrix=matrix)

            # If sorting by columns, transpose the matrix first
            if axis == "columns":
                matrix = self.transpose_matrix(matrix)

            # Separate the header from the data rows
            header = matrix[0]
            data_rows = matrix[1:]

            # Sort the data rows based on residue indices
            sorted_data_rows = sorted(data_rows, key=lambda row: int(row[0].replace(" ", "")[3:].split("-")[0]))

            # Combine the header with the sorted data rows
            sorted_matrix = [header] + sorted_data_rows

            # If sorting was by columns, transpose the sorted matrix back
            if axis == "columns":
                sorted_matrix = self.transpose_matrix(sorted_matrix)

            return sorted_matrix

        self._check_variable_types(
            variables=[
                interaction_data,
                axis,
                thr_interactions,
                thr_activity,
                selected_items,
                _count,
                residue_chain,
                save,
            ],
            expected_types=[
                InteractionData,
                str,
                (int, None.__class__),
                (float, None.__class__),
                (int, None.__class__),
                bool,
                bool,
                (str, None.__class__),
            ],
            variable_names=[
                "interaction_data",
                "axis",
                "thr_interactions",
                "thr_activity",
                "selected_items",
                "count",
                "residue_chain",
                "save",
            ],
        )

        data = copy.deepcopy(interaction_data)
        matrix = data.matrix

        self._verify_dimensions(matrix=matrix)

        if axis not in ["rows", "columns"]:
            raise InvalidAxisException(axis)

        # Raise an error if `thr_interactions` and `selected_items` are provided simultaneously
        if thr_interactions is not None and selected_items is not None:
            raise ValueError("Cannot select by both 'thr_interactions' and 'selected_items' at the same time.")
        if thr_interactions is not None and thr_activity is not None:
            raise ValueError("Cannot select by both 'thr_interactions' and 'thr_activity' at the same time.")
        if thr_activity is not None and selected_items is not None:
            raise ValueError("Cannot select by both 'thr_activity' and 'selected_items' at the same time.")

        # Ensure the correct axis is selected for activity-based selection
        if thr_activity is not None:
            axis = "rows" if self._get_residues_axis(matrix=matrix) == "columns" else "columns"

        # Transpose the matrix if operating on columns
        if axis == "columns":
            matrix = self.transpose_matrix(matrix)

        # Initialize a dictionary to store interaction counts per row/column
        reactives = {}
        for row in range(1, len(matrix)):
            for column in range(1, len(matrix[row])):
                cell = matrix[row][column]
                interactions = get_interactions(cell)
                reactives[row] = reactives.get(row, 0) + interactions

        # If `count` is True, return the interaction counts
        if _count:
            data = [list(reactives.keys()), list(reactives.values())]
            for index in reactives.keys():
                original_value = matrix[index][0]
                split_value = original_value.split("_")[0].strip()  # Splits and removes spaces
                data[0][index - 1] = split_value
            return data

        # Select rows/columns based on the provided criteria
        elif thr_interactions is not None:
            reactives = [
                key
                for key, value in sorted(reactives.items(), key=lambda item: item[1], reverse=True)
                if value >= thr_interactions
            ]
        elif self._get_residues_axis(matrix) == "columns" and thr_activity is not None:
            try:
                float(matrix[1][0].split(" (")[1].replace(")", ""))
            except Exception:
                raise MissedActivityException("The matrix does not contain activity values.")
            reactives = [
                key
                for key, value in sorted(
                    reactives.items(),
                    key=lambda item: float(matrix[item[0]][0].split(" (")[1].replace(")", "")),
                    reverse=True,
                )
                if float(matrix[key][0].split(" (")[1].replace(")", "")) >= thr_activity
            ]

        elif selected_items:
            selected_items = min(selected_items, len(matrix))
            reactives = [
                key for key, value in sorted(reactives.items(), key=lambda item: item[1], reverse=True)[:selected_items]
            ]
        else:
            reactives = [key for key, value in sorted(reactives.items(), key=lambda item: item[1], reverse=True)]

        # Create the selection matrix with the chosen rows/columns
        selection = [matrix[0]] + [matrix[row] for row in reactives]

        # Sort the selection by residue chain if specified
        if residue_chain:
            selection = sort_by_residue(matrix=matrix)

        # Transpose the selection back if it was initially transposed for columns
        if axis == "columns":
            selection = self.transpose_matrix(selection)

        data.matrix = selection

        if save:
            self.save_interaction_data(interaction_data=data, filename=save)

        return data

    def transpose_matrix(
        self, interaction_data: InteractionData, save: str = None, _pie: bool = False
    ) -> InteractionData:
        """
        Transposes the given interaction matrix.

        This method swaps rows and columns in the interaction matrix, effectively
        transposing it. If a save path is provided, the transposed matrix is stored
        as a file.

        Args:
            interaction_data (InteractionData): The interaction data containing the matrix to be transposed.
            save (str, optional): File path to save the transposed matrix. Defaults to None.

        Returns:
            InteractionData: The updated InteractionData object with the transposed matrix.

        Raises:
            TypeMismatchException: If the provided arguments have incorrect types.
            ValueError: If the matrix dimensions are invalid.
        """

        self._check_variable_types(
            variables=[interaction_data, save],
            expected_types=[(InteractionData, list), (str, None.__class__)],
            variable_names=["interaction_data", "save"],
        )

        data = copy.deepcopy(interaction_data)
        matrix = data.matrix if isinstance(interaction_data, InteractionData) else copy.deepcopy(data)

        self._verify_dimensions(matrix=matrix, _pie=_pie)

        # Transpose the matrix using list comprehension
        transposed = [[row[i] for row in matrix] for i in range(len(matrix[0]))]

        if isinstance(interaction_data, InteractionData):
            data.matrix = transposed
        else:
            data = transposed

        if save:
            self.save_interaction_data(interaction_data=data, filename=save)

        return data
