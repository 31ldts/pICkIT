"""Custom exceptions used across the pickit package.

Extracted verbatim from the original monolithic ``analyze_interactions.py``
(Fase 2, paso 1 del plan de modularización). No logic changes were made
during this extraction — only the physical location of the code changed.
"""


class TypeMismatchException(Exception):
    def __init__(self, variable_name, expected_types, actual_type):
        expected_types_str = ", ".join([t.__name__ for t in expected_types])
        self.message = (
            f"Variable '{variable_name}' has type {actual_type.__name__}, expected one of ({expected_types_str})."
        )
        super().__init__(self.message)


class FileOrDirectoryException(Exception):
    def __init__(self, path, error_type, message=None):
        self.path = path
        self.error_type = error_type  # Puede ser 'not_found' o 'empty'

        default_messages = {
            "not_found": f"File or directory not found: '{path}'",
            "empty": f"Directory '{path}' is empty",
        }

        self.message = message if message else default_messages.get(error_type, "An error occurred")
        super().__init__(self.message)


class InvalidColorException(Exception):
    """
    Exception raised when an invalid hexadecimal color is provided.
    """

    def __init__(self, invalid_colors: list[str]):
        self.invalid_colors = invalid_colors
        self.message = f"Invalid hexadecimal color(s) detected: {', '.join(invalid_colors)}"
        super().__init__(self.message)


class InvalidFileExtensionException(Exception):
    """Exception raised when the file extension is not .xlsx."""

    def __init__(self, filename, message="Invalid file extension. Only '.xlsx' files are allowed"):
        self.filename = filename
        self.message = f"{message}: '{filename}'"
        super().__init__(self.message)


class InvalidAxisException(Exception):
    """Exception raised for invalid axis values."""

    def __init__(self, axis_value):
        self.axis_value = axis_value
        self.message = f"Invalid axis value: '{axis_value}'. Expected 'rows' or 'columns'."
        super().__init__(self.message)


class InvalidFilenameException(Exception):
    """Exception raised for invalid mode values."""

    def __init__(self, filenames):
        self.filenames = filenames
        output = ""
        for filename in self.filenames:
            output += "\n\t- " + filename
        self.message = f"Some files have an invalid file name: {output}. \nFile names must not contain spaces."
        super().__init__(self.message)


class HeatmapActivityException(Exception):
    """Exception raised for invalid mode values."""

    def __init__(self):
        self.message = "Heatmap modes' max, min and mean require ligand/complex activities."
        super().__init__(self.message)


class InvalidModeException(Exception):
    """Exception raised for invalid mode values."""

    def __init__(self, mode, expected_values):
        self.message = f"Invalid mode value: '{mode}'. Expected:{expected_values}"
        super().__init__(self.message)


class MissedActivityException(Exception):
    """Exception raised when activity values are not found in the matrix."""

    def __init__(self, message):
        self.message = message
        super().__init__(self.message)
