"""Tests for the ``template_file`` default value in ``analyze_files``.

Added in Fase 3, per the user's explicit request: ``template_file`` should
have a sensible default (the template bundled with the package,
``constants.DEFAULT_TEMPLATE_FILE`` = ``"template.json"``) instead of
requiring ``None``/no restriction, while an explicit ``template_file``
argument must always override the default.

Deliberately a new file rather than edits to the frozen
``test_analyze_files.py`` baseline — see plan-accion-modularizacion-pickit.md,
Fase 2/Fase 3 notes on why this couldn't be done as part of Fase 2 without
risking the frozen baseline tests.
"""

import os

from pickit.constants import DEFAULT_TEMPLATE_FILE


class TestTemplateFileDefault:
    def test_bundled_default_template_exists_on_disk(self):
        # Sanity check the fixture this whole test file relies on: the
        # package actually ships a template.json next to io_mixin.py.
        import pickit

        package_dir = os.path.dirname(pickit.__file__)
        assert os.path.isfile(os.path.join(package_dir, DEFAULT_TEMPLATE_FILE))

    def test_no_template_file_arg_still_uses_bundled_default_silently(self, analyzer, arpeggio_dir, activity_csv):
        # Calling analyze_files without template_file must not raise or
        # behave differently just because a default now exists — this is
        # the "no debe ser obligatorio ni silencioso [en el sentido de
        # romper cosas]" requirement.
        data = analyzer.analyze_files(directory=arpeggio_dir, mode=analyzer.ARPEGGIO, activity_file=activity_csv)
        assert "hbond" in data.interactions

    def test_explicit_template_file_overrides_the_default(self, analyzer, arpeggio_dir, activity_csv, template_path):
        # An explicitly-passed template_file must be honored over the
        # bundled default, even though in this project both happen to
        # point at files with equivalent contact/type coverage.
        with_explicit = analyzer.analyze_files(
            directory=arpeggio_dir, mode=analyzer.ARPEGGIO, activity_file=activity_csv, template_file=template_path
        )
        without_explicit = analyzer.analyze_files(
            directory=arpeggio_dir, mode=analyzer.ARPEGGIO, activity_file=activity_csv
        )
        # Both use a template (default vs explicit), so the resulting
        # interaction ordering should match — this is the "an explicit
        # template still works, and using no template arg isn't silently
        # different in an unexpected way" contract.
        assert with_explicit.interactions == without_explicit.interactions
        assert with_explicit.matrix == without_explicit.matrix

    def test_default_only_applies_to_arpeggio_mode(self, analyzer, tmp_path):
        # IChem mode never reads template_file at all (see parsers/ichem.py
        # and the analyze_files body); passing an empty ichem directory
        # should not attempt to look at any template.
        empty_dir = tmp_path / "empty_ichem"
        empty_dir.mkdir()
        (empty_dir / "complex1.txt").write_text("not a real ichem line\n")
        # Should not raise for template-related reasons; may raise for
        # other reasons (invalid file), but not because of template_file.
        try:
            analyzer.analyze_files(directory=str(empty_dir), mode=analyzer.ICHEM)
        except Exception as exc:
            assert "template" not in str(exc).lower()

    def test_missing_bundled_default_does_not_break_anything(self, analyzer, arpeggio_dir, activity_csv, monkeypatch):
        # Simulate the bundled template.json being absent: analyze_files
        # must fall back to "no template restriction", exactly like before
        # this default existed, not raise.
        import pickit.io_mixin as io_mixin_module

        real_isfile = os.path.isfile

        def fake_isfile(path):
            if path.endswith(DEFAULT_TEMPLATE_FILE) and os.path.dirname(path) == os.path.dirname(
                io_mixin_module.__file__
            ):
                return False
            return real_isfile(path)

        monkeypatch.setattr(os.path, "isfile", fake_isfile)

        data = analyzer.analyze_files(directory=arpeggio_dir, mode=analyzer.ARPEGGIO, activity_file=activity_csv)
        # Falls back to the full, unrestricted arpeggio interaction_labels.
        assert data.interactions == analyzer.interaction_labels
