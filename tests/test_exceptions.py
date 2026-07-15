"""Tests for pickit.exceptions.

New tests introduced during Fase 2, step 1 (extraction of exceptions.py).
The exception *behavior* (when each is raised) is already covered by the
Fase 0 baseline suite through the public API (test_configuration.py,
test_analyze_files.py, etc.) — these tests only check the exceptions'
own construction/message logic in isolation, which had no dedicated
coverage before.
"""

import pytest

from pickit import exceptions


def test_type_mismatch_exception_message():
    exc = exceptions.TypeMismatchException("plot_max_cols", (int,), str)
    assert "plot_max_cols" in exc.message
    assert "int" in exc.message
    assert "str" in exc.message


def test_file_or_directory_exception_default_messages():
    not_found = exceptions.FileOrDirectoryException("foo", "not_found")
    assert "not found" in not_found.message
    empty = exceptions.FileOrDirectoryException("foo", "empty")
    assert "empty" in empty.message


def test_invalid_color_exception_lists_colors():
    exc = exceptions.InvalidColorException(["not-a-color", "also-bad"])
    assert "not-a-color" in exc.message
    assert "also-bad" in exc.message


def test_invalid_file_extension_exception():
    exc = exceptions.InvalidFileExtensionException("out.csv")
    assert "out.csv" in exc.message


def test_invalid_axis_exception():
    exc = exceptions.InvalidAxisException("diagonal")
    assert "diagonal" in exc.message


def test_invalid_filename_exception_lists_all_filenames():
    exc = exceptions.InvalidFilenameException(["bad name.json", "other bad.json"])
    assert "bad name.json" in exc.message
    assert "other bad.json" in exc.message


def test_heatmap_activity_exception_static_message():
    exc = exceptions.HeatmapActivityException()
    assert "activities" in exc.message


def test_invalid_mode_exception():
    exc = exceptions.InvalidModeException("sideways", ["input", "output"])
    assert "sideways" in exc.message


def test_missed_activity_exception_passthrough_message():
    exc = exceptions.MissedActivityException("custom message")
    assert exc.message == "custom message"


@pytest.mark.parametrize(
    "name",
    [
        "TypeMismatchException",
        "FileOrDirectoryException",
        "InvalidColorException",
        "InvalidFileExtensionException",
        "InvalidAxisException",
        "InvalidFilenameException",
        "HeatmapActivityException",
        "InvalidModeException",
        "MissedActivityException",
    ],
)
def test_reexported_from_analyze_interactions_is_same_class(name):
    # Backward compatibility: `from pickit.analyze_interactions import X`
    # (as tests/conftest.py already does) must resolve to the exact same
    # class object as `from pickit.exceptions import X`.
    import pickit.analyze_interactions as ai

    assert getattr(ai, name) is getattr(exceptions, name)
