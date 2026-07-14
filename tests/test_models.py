"""Tests for pickit.models.InteractionData.

New tests introduced during Fase 2, step 2 (extraction of models.py).
InteractionData usage as a data carrier is already exercised indirectly by
almost every test in the Fase 0 baseline suite; these tests focus on the
one piece of actual logic it has: `compare`.
"""

import copy

from pickit.models import InteractionData


def _make(**overrides):
    defaults = dict(
        colors=["#e6194B"],
        interactions=["hbond"],
        ligand=True,
        matrix=[["", "LIGA"], ["HIS 41", ""]],
        mode="arpeggio",
        protein=True,
        subunit=False,
        subunits_set={"A"},
    )
    defaults.update(overrides)
    return InteractionData(**defaults)


def test_compare_identical_objects_reports_no_changes():
    a = _make()
    b = copy.deepcopy(a)
    assert a.compare(b) == "There are no changes."


def test_compare_reports_each_differing_attribute():
    a = _make()
    b = _make(mode="ichem", protein=False)
    diff = a.compare(b)
    assert diff == {
        "mode": ("arpeggio", "ichem"),
        "protein": (True, False),
    }


def test_compare_against_object_without_dict_is_invalid():
    a = _make()
    assert a.compare(42) == "Invalid comparison. The object does not has the espected attributes."


def test_constructor_stores_all_fields_as_given():
    a = _make()
    assert a.colors == ["#e6194B"]
    assert a.interactions == ["hbond"]
    assert a.ligand is True
    assert a.matrix == [["", "LIGA"], ["HIS 41", ""]]
    assert a.mode == "arpeggio"
    assert a.protein is True
    assert a.subunit is False
    assert a.subunits_set == {"A"}


def test_reexported_from_analyze_interactions_is_same_class():
    # Backward compatibility: `from pickit.analyze_interactions import
    # InteractionData` (as tests/conftest.py and pickit/__init__.py do)
    # must resolve to the exact same class object.
    from pickit import analyze_interactions as ai

    assert ai.InteractionData is InteractionData
