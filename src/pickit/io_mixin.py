"""Directory/config management and file parsing (Arpeggio, IChem).

Extracted verbatim from the original monolithic ``analyze_interactions.py``
(Fase 2, paso 7 del plan de modularización — el de mayor riesgo). No logic
changes were made during this extraction — only the physical location of
the code changed.

Deliberately NOT split further into ``parsers/arpeggio.py`` +
``parsers/ichem.py`` in this step, even though that's what Anexo A
ultimately maps them to. ``ichem_analysis``, ``arpeggio_analysis`` and
``arpeggio_analysis_template`` are nested closures inside ``analyze_files``
that capture ``protein``/``ligand``/``subunit`` from the enclosing call and
several module-level constants; pulling them out into standalone modules
means turning that implicit closure state into explicit parameters, which
is a real refactor, not a "copiar-pegar primero" move. Per Anexo B's own
practical criterion ("no cambiar la estructura solo por modernizar" / no
unjustified redesign), that split is left as an explicit follow-up
sub-task rather than done silently here, alongside the ``template_file``
default-value sub-task (see plan-accion-modularizacion-pickit.md) — both
touch this same code and are best done together, with their own dedicated
test coverage, rather than folded into a same-behavior extraction step.
"""

import csv
import json
import os
import re

from .constants import (
    AMINO_ACID_CODES,
    ARPEGGIO_COLORS,
    ARPEGGIO_CONT,
    ARPEGGIO_INT_ENT,
    ARPEGGIO_TYPE,
    COLORS,
    DIFF_DELIM,
    GROUP_DELIM,
    INTERACTION_LABELS,
    PROGRAM_MODES,
    SAME_DELIM,
)
from .exceptions import FileOrDirectoryException, InvalidColorException, InvalidModeException
from .models import InteractionData


class IOMixin:
    """Directory configuration and interaction-file parsing.

    Depends on ``ValidationMixin`` (``self._check_variable_types``),
    ``FilterMixin`` (``self.sort_matrix``), ``ExportMixin``
    (``self.save_interaction_data``), and on configuration state
    (``self.input_directory``, ``self.saving_directory``,
    ``self.interaction_labels``, ``self.plot_colors``, ``self.aa``).
    """

    def _logger(self, type: str, message: str) -> None:
        """
        Displays a message to the user.

        Args:
            type (str): Type of message ('info', 'warning', 'error').
            message (str): The message to display.

        Returns:
            None
        """
        if type == "info":
            print(f"\033[1;34mINFO\033[0m: {message}")  # Azul (34)
        elif type == "warning":
            print(f"\033[1;33mWARNING\033[0m: {message}")  # Amarillo (33)
        elif type == "error":
            print(f"\033[1;31mERROR\033[0m: {message}")  # Rojo (31)
        else:
            raise ValueError("Invalid message type. Expected 'info', 'warning', or 'error'.")

    def change_directory(self, path: str, mode: str) -> None:
        """
        Changes the working directory for saving or input operations.

        Args:
            path (str): Name of the subdirectory to switch to.
            mode (str): Determines whether to change the input or output directory.
                - 'input': Sets the directory for input files.
                - 'output': Sets the directory for output files.

        Returns:
            None

        Raises:
            ValueError: If the specified directory does not exist.
            InvalidModeException: If an invalid mode is provided.
        """

        self._check_variable_types(variables=[path, mode], expected_types=[str, str], variable_names=["path", "mode"])

        # Construct the new path
        new_path = os.path.join(os.getcwd(), path)

        # Verify the new path exists
        if not os.path.exists(new_path):
            raise ValueError("The specified directory must exist inside the project.")

        if mode == "input":
            self.input_directory = new_path
        elif mode == "output":
            self.saving_directory = new_path
        else:
            raise InvalidModeException(mode=mode, expected_values=["input", "output"])

    def set_config(
        self,
        interactions: list[str] = None,
        plot_max_cols: int = None,
        plot_colors: list[str] = None,
        reset: bool = False,
        mode: str = None,
        heat_max_cols: int = None,
        heat_colors: str = None,
        interaction_data: InteractionData = None,
    ) -> None:
        """
        Configures interaction settings, including labels, colors, and visualization parameters.

        Args:
            interactions (list[str], optional): List of interaction labels.
            plot_max_cols (int, optional): Maximum number of columns for plot visualization.
            plot_colors (list[str], optional): List of colors in hexadecimal format.
            reset (bool, optional): If True, resets configurations to default values.
            mode (str, optional): Determines preset configurations for different analysis modes ('ichem' or 'arpeggio').
            heat_max_cols (int, optional): Maximum number of columns for heatmap visualization.
            heat_colors (str, optional): Color scheme for heatmaps.
            interaction_data (InteractionData, optional): Object containing interaction settings to be applied.

        Returns:
            None

        Raises:
            InvalidColorException: If any color in the provided list is not a valid hexadecimal value.
            InvalidModeException: If an invalid mode is provided.
        """

        def is_valid_hex_color(color: str) -> bool:
            """
            Validates if a string is a valid hexadecimal color.

            Args:
                color (str): Color string to validate.

            Returns:
                bool: True if the color is a valid hex value, False otherwise.
            """
            return bool(re.match(r"^#([0-9a-fA-F]{6}|[0-9a-fA-F]{3})$", color))

        def reset_configuration() -> None:
            """
            Resets interaction labels and colors to their default values.

            Returns:
                None
            """
            self.saving_directory = os.getcwd()  # Set the default saving directory
            self.input_directory = os.getcwd()  # Set the default input directory
            self.interaction_labels = INTERACTION_LABELS  # Default interaction labels
            self.codes = True
            self.plot_colors = COLORS  # Default color configuration
            self.plot_max_cols = 80
            self.aa = AMINO_ACID_CODES
            self.heat_max_cols = 30
            self.heat_colors = "RdYlGn"

        self._check_variable_types(
            variables=[
                interactions,
                plot_max_cols,
                plot_colors,
                reset,
                mode,
                heat_max_cols,
                heat_colors,
                interaction_data,
            ],
            expected_types=[
                (list, None.__class__),
                (int, None.__class__),
                (list, None.__class__),
                bool,
                (str, None.__class__),
                (int, None.__class__),
                (str, None.__class__),
                (InteractionData, None.__class__),
            ],
            variable_names=[
                "interactions",
                "plot_max_cols",
                "colors",
                "reset",
                "mode",
                "heat_max_elems",
                "heat_colors",
                "interaction_data",
            ],
        )

        # Perform reset if requested
        if reset:
            reset_configuration()
            return

        # Update interaction labels if provided
        if interactions:
            self.interaction_labels = interactions

        # Update colors if provided and valid
        if plot_colors:
            invalid_colors = [color for color in plot_colors if not is_valid_hex_color(color)]
            if invalid_colors:
                raise InvalidColorException(invalid_colors)
            self.plot_colors = plot_colors

        if mode:
            if mode == "ichem":
                self.interaction_labels = INTERACTION_LABELS  # Default interaction labels
                self.plot_colors = COLORS  # Default color configuration
            elif mode == "arpeggio":
                self.interaction_labels = ARPEGGIO_CONT[:2] + ARPEGGIO_TYPE + ARPEGGIO_CONT[2:]
                # self.interaction_labels.sort(key=lambda s: (s.lower(), s.islower()))
                self.plot_colors = ARPEGGIO_COLORS  # Default color c
            else:
                raise InvalidModeException(mode=mode, expected_values=PROGRAM_MODES)
        if plot_max_cols:
            self.plot_max_cols = plot_max_cols
        if heat_max_cols:
            self.heat_max_cols = heat_max_cols
        if heat_colors:
            self.heat_colors = heat_colors
        if interaction_data:
            self.interaction_labels = interaction_data.interactions
            self.plot_colors = interaction_data.colors

    def analyze_files(
        self,
        directory: str,
        mode: str,
        activity_file: str = None,
        protein: bool = True,
        ligand: bool = True,
        subunit: bool = False,
        template_file: str = None,
        save: str = None,
    ) -> InteractionData:
        """
        Analyzes interaction data files in a specified directory, processing them according to the specified mode.

        This function reads interaction data files, extracts relevant interaction details, and organizes them into a
        structured matrix format. Optionally, an activity file can be provided to label data based on activity levels.

        Args:
            directory (str): Path to the directory containing interaction data files.
            mode (str): Processing mode. Supported modes:
                - 'ichem': Processes IChem interaction files.
                - 'arpeggio': Processes Arpeggio interaction files.
            activity_file (str, optional): Path to a CSV file containing activity data for annotation.
            protein (bool, optional): Whether to include protein atoms in the analysis. Defaults to True.
            ligand (bool, optional): Whether to include ligand atoms in the analysis. Defaults to True.
            subunit (bool, optional): Whether to differentiate between protein subunits. Defaults to False.
            template_file (str, optional): Path to a JSON file containing a template for the interaction data.
            save (str, optional): Path to save the processed interaction matrix. Defaults to None.

        Returns:
            InteractionData: An object containing the processed interaction matrix, metadata, and interaction labels.

        Raises:
            FileNotFoundError: If the specified directory or activity file does not exist.
            EmptyDirectoryException: If the specified directory contains no valid files.
            InvalidModeException: If the provided mode is not supported.
            InvalidFilenameException: If any file names contain spaces or invalid characters.
        """

        def label_matrix(
            matrix: list[list[str]],
            rows: list[str],
            columns: list[str],
            activity_file: str,
            correction: list[str] = None,
        ) -> list[list[str]]:
            """
            Adds appropriate headers to the interaction matrix, including residue and file names.

            If an activity file is provided, it also labels the columns with corresponding activity values.

            Args:
                matrix (list[list[str]]): 2D list representing interaction data.
                rows (list[str]): List of residue names for row labeling.
                columns (list[str]): List of file names for column labeling.
                activity_file (str): Path to the activity file for activity-based labeling.
                correction (list[str], optional): Optional list to correct column names.

            Returns:
                list[list[str]]: The labeled interaction matrix.
            """
            rows = [row.replace("\t", "") for row in rows]
            if self.codes and correction:
                columns = [correction[cont] + ", " + column for cont, column in enumerate(columns)]
            columns = [""] + columns[:]  # Add an empty string at the start for residue names

            if activity_file:
                if not os.path.isfile(activity_file):
                    raise FileNotFoundError(f"The file '{activity_file}' does not exist.")

                # Read activity data into a dictionary
                data_dict = {}
                with open(activity_file, newline="") as csvfile:
                    csvreader = csv.reader(csvfile)
                    try:
                        next(csvreader)  # Skip header
                    except Exception:
                        raise ValueError(f"The CSV file '{activity_file}' is missing a header.")
                    for key, value in csvreader:
                        data_dict[key.upper()] = str(round(float(value), 3))

                if not data_dict:
                    raise ValueError(f"The CSV file '{activity_file}' must contain at least one row of data.")

                # Update column names with activity data
                unprocessed_files = 0
                if correction:
                    for i in range(1, len(columns)):
                        drug_name = columns[i]
                        activity = data_dict.get(correction[i - 1], "0")
                        unprocessed_files += 1 if activity == "0" else 0
                        columns[i] = f"{drug_name} ({activity})"
                else:
                    for i in range(1, len(columns)):
                        drug_name = columns[i]
                        activity = data_dict.get(drug_name, "0")
                        unprocessed_files += 1 if activity == "0" else 0
                        columns[i] = f"{drug_name} ({activity})"

                if unprocessed_files > 0:
                    self._logger(
                        type=self.WARNING,
                        message=(
                            f"{unprocessed_files} ({unprocessed_files * 100 / (len(columns) - 1):.2f} %) "
                            "files were not in the activity file."
                        ),
                    )
                    self._logger(
                        type=self.INFO,
                        message=(
                            "Files not found in the activity file are labeled as '0'. "
                            "Use sort_matrix(interaction_data=data, thr_activity=0.001) to filter them out."
                        ),
                    )
                else:
                    self._logger(type=self.INFO, message="All files were in the activity file.")

            else:
                for i in range(1, len(columns)):
                    drug_name = columns[i]
                    columns[i] = f"{drug_name}"

            # Insert headers into the matrix
            matrix.insert(0, columns)
            for i, row in enumerate(matrix[1:], start=1):
                row.insert(0, rows[i - 1])

            return matrix

        def modify_cell(
            text: str,
            interaction: str,
            atoms: str,
            interaction_labels: list,
        ) -> str:
            """
            Updates the cell content by adding interaction type and involved atoms.

            Args:
                text (str): Current content of the matrix cell.
                interaction (str): Type of interaction occurring.
                atoms (str): Atoms involved in the interaction.
                interaction_labels (list): List of predefined interaction labels.

            Returns:
                str: Updated cell content with formatted interaction information.
            """
            # Create an interaction_map based on the global INTERACTION_LABELS list
            interaction_map = {label: str(index + 1) for index, label in enumerate(interaction_labels)}

            # Assign the interaction code based on the interaction_map, or default to the last value
            interaction_code = interaction_map.get(interaction, str(len(interaction_labels)))

            # If the text is empty, add the interaction with the provided atoms
            if text == "":
                return f"{interaction_code} {GROUP_DELIM}{atoms}{GROUP_DELIM}"

            # Split the cell content and remove empty parts, keeping existing interactions
            content = text.replace(DIFF_DELIM, "").split(GROUP_DELIM)[:-1]
            exists = False
            cell = ""

            # Check if the interaction already exists and add atoms to the corresponding interaction
            for index, segment in enumerate(content):
                if index % 2 == 0 and interaction_code == segment.strip():
                    # Delete repate entries
                    entries = content[index + 1].split(SAME_DELIM)
                    detected = False
                    for entry in entries:
                        if entry == atoms:
                            detected = True
                            break
                    if not detected:
                        content[index + 1] += f"{SAME_DELIM}{atoms}"
                    exists = True
                    break

            # Rebuild the cell content, preserving existing interactions
            cell = DIFF_DELIM.join(
                f"{content[i]}{GROUP_DELIM}{content[i + 1]}{GROUP_DELIM}" for i in range(0, len(content), 2)
            )

            # If the interaction didn't exist, append it at the end
            if not exists:
                cell += f"{DIFF_DELIM}{interaction_code} {GROUP_DELIM}{atoms}{GROUP_DELIM}"

            return cell

        def read_file(file_name: str) -> list[str]:
            """
            Reads the specified file and returns its content as a list of lines.

            Supports reading both JSON and text-based interaction files.

            Args:
                file_name (str): The name of the file to read.

            Returns:
                list[str]: List of lines from the file (for text files) or parsed JSON data.

            Raises:
                FileNotFoundError: If the file does not exist.
                Exception: If an unexpected error occurs during reading.
            """
            try:
                with open(file_name, "r") as file:
                    if file_name.split(".")[-1] == "json":
                        return json.load(file)
                    else:
                        return [line.strip() for line in file.readlines()]
            except FileNotFoundError:
                print(f"Error: The file '{file_name}' does not exist.")
                raise
            except Exception as e:
                print(f"Error: An unexpected error occurred while reading '{file_name}': {e}")
                raise

        def adjust_subunits(
            matrix: list[list[str]],
        ) -> list[list[str]]:
            """
            Adjusts matrix entries to correctly reflect subunit information for residues.

            Removes duplicate atoms and ensures that subunit data is accurately represented.

            Args:
                matrix (list[list[str]]): The interaction matrix to modify.

            Returns:
                list[list[str]]: Updated matrix with subunit adjustments.
            """

            def remove_duplicate_atoms(atoms: str) -> str:
                """
                Removes duplicate atoms from the input string.

                Args:
                    atoms (str): A string containing atom interactions separated by SAME_DELIM.

                Returns:
                    str: A formatted string with unique atoms, enclosed in '|' delimiters.
                """
                # Split the input string by SAME_DELIM to get individual atoms.
                sections = atoms.split(SAME_DELIM)

                # Use a set to keep only unique atoms in their original order.
                unique_atoms = list(dict.fromkeys(sections))

                # Join the unique atoms back into a single string separated by SAME_DELIM.
                text = SAME_DELIM.join(unique_atoms)

                # Return the formatted string enclosed in '|' delimiters.
                return f"{GROUP_DELIM}{text}{GROUP_DELIM}"

            for row in range(len(matrix)):
                for column in range(len(matrix[row])):
                    cell = matrix[row][column]
                    if cell != "":
                        sections = cell.split(GROUP_DELIM)
                        text = "".join(
                            sections[i - 1] + remove_duplicate_atoms(sections[i]) for i in range(1, len(sections), 2)
                        )
                        matrix[row][column] = text
            return matrix

        def sort_interactions(matrix: list[list[str]]) -> list[list[str]]:
            """
            Sorts the interaction types within each cell of the matrix in ascending order.

            Args:
                matrix (list[list[str]]): The matrix to be sorted.

            Returns:
                list[list[str]]: The sorted matrix with interactions ordered numerically.
            """
            for row_index, row in enumerate(matrix):
                for cell_index, cell in enumerate(row):
                    if cell != "":
                        # Split the cell into individual interactions
                        interactions = cell.split(DIFF_DELIM)

                        # Sort the interactions based on the number that follows the initial space
                        if len(interactions) > 1:
                            interactions = sorted(interactions, key=lambda x: int(x.split(" ")[0]))

                        # Join the sorted interactions back and update the cell
                        matrix[row_index][cell_index] = DIFF_DELIM.join(interactions)

            return matrix

        def validate_string(input_string: str) -> bool:
            """
            Validates a residue string to ensure it follows the expected format.

            Format:
            - Starts with a three-letter amino acid abbreviation.
            - Followed by spaces and a numeric sequence.
            - Ends with a dash and additional characters.

            Args:
                input_string (str): The residue string to validate.

            Returns:
                bool: True if the format is valid, False otherwise.
            """
            # Regular expression to validate the required format
            pattern = r"^[A-Z]{3} +\d+-.+$"

            # Match the input string with the regular expression pattern
            return bool(re.match(pattern, input_string))

        def get_protein_ligand(begin: dict, end: dict) -> tuple[dict, dict]:
            if begin["label_comp_type"] == "P" and end["label_comp_type"] == "P":
                if begin["label_comp_id"] in self.aa and end["label_comp_id"] in self.aa:
                    return None, None
                elif begin["label_comp_id"] in self.aa:
                    return begin, end
                else:
                    return end, begin
            elif begin["label_comp_type"] == "P":
                return begin, end
            elif end["label_comp_type"] == "P":
                return end, begin
            else:
                return None, None

        def validate_file(filename, mode):
            """
            Checks if a filename is valid (i.e., contains no spaces).

            Args:
                filename (str): Name of the file.

            Returns:
                bool: True if the filename is valid, False otherwise.
            """
            if filename.count(" ") != 0:
                print(f"Warning: The filename '{filename}' contains spaces.")
                return False
            if mode == self.ICHEM:
                if filename.split(".")[-1] != "txt":
                    print(f"Warning: The filename '{filename}' is not a valid IChem file.")
                    return False
            elif mode == self.ARPEGGIO:
                if filename.split(".")[-1] != "json":
                    print(f"Warning: The filename '{filename}' is not a valid Arpeggio file.")
                    return False
            return True

        def check_directory(directory):
            """
            Validates and retrieves the list of files in the specified directory.

            Args:
                directory (str): Path to the directory.

            Returns:
                list[str]: List of filenames in the directory.

            Raises:
                FileNotFoundError: If the directory does not exist.
                EmptyDirectoryException: If the directory is empty.
            """
            if not os.path.exists(directory):
                raise FileOrDirectoryException(path=directory, error_type="not_found")
            elif not os.path.isdir(directory):
                raise FileOrDirectoryException(path=directory, error_type="not_found")
            else:
                files = os.listdir(directory)
                if not files:
                    raise FileOrDirectoryException(path=directory, error_type="empty")
            return files

        def ichem_analysis(content, index, files, subunits_set, cont, matrix, aa):
            """
            Processes IChem interaction files, extracting relevant data and updating the matrix.
            """
            for line in content:
                elements = line.split(GROUP_DELIM)
                if len(elements) == 10:
                    interaction = elements[0].strip().replace("\t", "")
                    residue = elements[3].strip().replace("\t", "")
                    if validate_string(residue):
                        if not subunit:
                            sections = residue.split("-")
                            residue = sections[0]
                            subunits_set.add(sections[1])

                        atoms = (
                            f"{elements[1].strip()}-{elements[4].strip()}"
                            if protein and ligand
                            else elements[1].strip()
                            if protein
                            else elements[4].strip()
                            if ligand
                            else ""
                        )
                        if not subunit:
                            atoms += f"({sections[1]})"

                        if residue not in aa:
                            aa[residue] = cont
                            cont += 1
                        column = aa[residue]

                        # Ensure matrix size and modify cell
                        if len(matrix) <= column:
                            matrix.append([""] * len(files))

                        matrix[column][index] = modify_cell(
                            text=matrix[column][index],
                            interaction=interaction,
                            atoms=atoms,
                            interaction_labels=INTERACTION_LABELS,
                        )
            return matrix, aa, cont, subunits_set

        def arpeggio_analysis(content, index, files, subunits_set, cont, matrix, aa):
            """
            Processes Arpeggio interaction files, extracting relevant data and updating the matrix.
            """
            interaction_list = ARPEGGIO_CONT[:2] + ARPEGGIO_TYPE + ARPEGGIO_CONT[2:]
            # interaction_list.sort(key=lambda s: (s.lower(), s.islower()))

            # Filter to obtain entries with interacting_entities == INTER
            inter_set = [elem for elem in content if elem["interacting_entities"] in ARPEGGIO_INT_ENT]
            # Filter to obtain entries with the desired contact or type
            inter_set = [
                elem
                for elem in inter_set
                for contact in elem["contact"]
                if contact in ARPEGGIO_CONT or elem["type"] in ARPEGGIO_TYPE
            ]
            for inter in inter_set:
                if inter["type"] in ARPEGGIO_TYPE:
                    contact = [inter["type"]]
                else:
                    contact = [conta for conta in inter["contact"] if conta in ARPEGGIO_CONT]
                prot, lig = get_protein_ligand(begin=inter["bgn"], end=inter["end"])
                if prot is not None:
                    residue = prot["label_comp_id"] + " " + str(prot["auth_seq_id"])
                    prot_atom = prot["auth_atom_id"]
                    prot_subunit = prot["auth_asym_id"]
                    ligand_code = lig["label_comp_id"]
                    lig_atom = lig["auth_atom_id"]

                    subunits_set.add(prot_subunit)
                    atoms = (
                        f"{prot_atom}-{lig_atom}"
                        if protein and ligand
                        else prot_atom
                        if protein
                        else lig_atom
                        if ligand
                        else ""
                    )

                    if subunit:
                        residue += "-" + prot_subunit
                    else:
                        atoms += f"({prot_subunit})"

                    if residue not in aa:
                        aa[residue] = cont
                        cont += 1
                    column = aa[residue]

                    # Ensure matrix size and modify cell
                    if len(matrix) <= column:
                        matrix.append([""] * len(files))

                    for interaction in contact:
                        matrix[column][index] = modify_cell(
                            text=matrix[column][index],
                            interaction=interaction,
                            atoms=atoms,
                            interaction_labels=interaction_list,
                        )

            return matrix, ligand_code, aa, cont, subunits_set

        def arpeggio_analysis_template(
            content, index, files, subunits_set, cont, matrix, aa, template, interaction_list
        ):
            """
            Processes Arpeggio interaction files, extracting relevant data and updating the matrix.
            """

            def matches_template(entry, template):
                """
                Compara una entrada con un template recursivamente.
                """
                if isinstance(template, dict):
                    if not isinstance(entry, dict):
                        return False
                    for key, tmpl_val in template.items():
                        if key not in entry:
                            return False
                        if not matches_template(entry[key], tmpl_val):
                            return False
                    return True
                elif isinstance(template, list):
                    if template is None:
                        return True
                    # Aquí asumimos que entry debe ser un valor único que esté en la lista
                    return any(value in template for value in entry)
                elif template is None:
                    return True
                else:
                    return entry == template

            # Filter to obtain entries with the desired contact or type
            for inter in content:
                for entry in template:
                    if matches_template(inter, entry):
                        contact = set()
                        if type(inter["type"]) == list:
                            for conta in inter["type"]:
                                if conta in interaction_list:
                                    contact.add(conta)
                        else:
                            if inter["type"] in interaction_list:
                                contact.add(inter["type"])
                        if type(inter["contact"]) == list:
                            for conta in inter["contact"]:
                                if conta in interaction_list:
                                    contact.add(conta)
                        else:
                            if inter["contact"] in interaction_list:
                                contact.add(inter["contact"])
                        prot, lig = get_protein_ligand(begin=inter["bgn"], end=inter["end"])
                        if prot is not None:
                            residue = prot["label_comp_id"] + " " + str(prot["auth_seq_id"])
                            prot_atom = prot["auth_atom_id"]
                            prot_subunit = prot["auth_asym_id"]
                            ligand_code = lig["label_comp_id"]
                            lig_atom = lig["auth_atom_id"]

                            subunits_set.add(prot_subunit)
                            atoms = (
                                f"{prot_atom}-{lig_atom}"
                                if protein and ligand
                                else prot_atom
                                if protein
                                else lig_atom
                                if ligand
                                else ""
                            )

                            if subunit:
                                residue += "-" + prot_subunit
                            else:
                                atoms += f"({prot_subunit})"

                            if residue not in aa:
                                aa[residue] = cont
                                cont += 1
                            column = aa[residue]

                            # Ensure matrix size and modify cell
                            if len(matrix) <= column:
                                matrix.append([""] * len(files))

                            for interaction in contact:
                                matrix[column][index] = modify_cell(
                                    text=matrix[column][index],
                                    interaction=interaction,
                                    atoms=atoms,
                                    interaction_labels=interaction_list,
                                )
                        break

            return matrix, ligand_code, aa, cont, subunits_set

        # Validate input types
        self._check_variable_types(
            variables=[directory, mode, activity_file, protein, ligand, subunit, template_file, save],
            expected_types=[
                str,
                str,
                (str, None.__class__),
                bool,
                bool,
                bool,
                (str, None.__class__),
                (str, None.__class__),
            ],
            variable_names=[
                "directory",
                "mode",
                "activity_file",
                "protein",
                "ligand",
                "subunit",
                "template_file",
                "save",
            ],
        )

        directory = os.path.join(self.input_directory, directory)
        if activity_file is not None:
            activity_file = os.path.join(self.input_directory, activity_file)
        if template_file is not None:
            template_file = os.path.join(self.input_directory, template_file)

        # Check the directory and return its files
        files = check_directory(directory=directory)

        # Check if the mode is registered
        if mode not in PROGRAM_MODES:
            raise InvalidModeException(mode=mode, expected_values=PROGRAM_MODES)
        else:
            self.set_config(mode=mode)

        ligands = [None] * len(files)
        matrix = []
        aa = {}
        cont = 0
        subunits_set = set()
        failed_files = []

        interaction_list = None
        if template_file is not None:
            with open(template_file) as f:
                template = json.load(f)
            # Get set of interactions from the template
            interaction_set = set()
            for entry in template:
                for field in ("contact", "type"):
                    if field in entry and entry[field] is not None:
                        if isinstance(entry[field], str):
                            interaction_set.add(entry[field])
                        else:
                            for interaction in entry[field]:
                                interaction_set.add(interaction)
            interaction_list = list(interaction_set)
            global_order = ARPEGGIO_CONT[:2] + ARPEGGIO_TYPE + ARPEGGIO_CONT[2:]
            # Filtra solo los que están presentes en el set y mantén el orden de global_order
            interaction_list = [item for item in global_order if item in interaction_set]
            # interaction_list.sort(key=lambda s: (s.lower(), s.islower()))

        # Analyze each file in the directory
        for index, file in enumerate(files):
            file_path = os.path.join(directory, file)

            if os.path.isfile(file_path) and validate_file(file, mode):
                content = read_file(file_path)
                if mode == self.ICHEM:
                    matrix, aa, cont, subunits_set = ichem_analysis(
                        content=content,
                        index=index,
                        files=files,
                        subunits_set=subunits_set,
                        cont=cont,
                        matrix=matrix,
                        aa=aa,
                    )
                    ligands[index] = file.replace(".txt", "").upper()
                elif mode == self.ARPEGGIO:
                    if template_file is not None:
                        matrix, ligand_code, aa, cont, subunits_set = arpeggio_analysis_template(
                            content=content,
                            index=index,
                            files=files,
                            subunits_set=subunits_set,
                            cont=cont,
                            matrix=matrix,
                            aa=aa,
                            template=template,
                            interaction_list=interaction_list,
                        )
                    else:
                        matrix, ligand_code, aa, cont, subunits_set = arpeggio_analysis(
                            content=content,
                            index=index,
                            files=files,
                            subunits_set=subunits_set,
                            cont=cont,
                            matrix=matrix,
                            aa=aa,
                        )
                    files[index] = file.replace(".json", "").upper()
                    ligands[index] = ligand_code

            else:
                failed_files.append(file_path)

        if mode == self.ICHEM:
            files = None

        if len(failed_files) > 0:
            files = [f for f in files if f not in failed_files]
            ligands = [ligand_ for ligand_ in ligands if ligand_ is not None]

        if not subunit:
            matrix = adjust_subunits(matrix=matrix)
        matrix = sort_interactions(matrix=matrix)
        matrix = label_matrix(
            matrix=matrix, rows=list(aa.keys()), columns=ligands, activity_file=activity_file, correction=files
        )

        interactions = interaction_list if interaction_list is not None else self.interaction_labels
        if len(self.plot_colors) > len(interactions):
            colors = self.plot_colors[: len(interactions)]
        elif len(self.plot_colors) < len(interactions):
            colors = self.plot_colors + [self.plot_colors[-1]] * (len(interactions) - len(self.plot_colors))
        else:
            colors = self.plot_colors.copy()

        interaction_data = InteractionData(
            colors=colors,
            interactions=interactions,
            ligand=ligand,
            matrix=matrix,
            mode=mode,
            protein=protein,
            subunit=subunit,
            subunits_set=subunits_set,
        )
        interaction_data = self.sort_matrix(interaction_data=interaction_data, residue_chain=True)
        # Save the matrix if specified
        if save:
            self.save_interaction_data(interaction_data=interaction_data, filename=save)

        return interaction_data
