import os

import pytest

from pickit.analyze_interactions import FileOrDirectoryException, InvalidModeException


class TestAnalyzeFilesArpeggio:
    def test_happy_path_builds_expected_matrix(self, analyzer, arpeggio_dir, activity_csv):
        data = analyzer.analyze_files(
            directory=arpeggio_dir,
            mode=analyzer.ARPEGGIO,
            activity_file=activity_csv,
        )

        header = data.matrix[0]
        row_labels = [row[0] for row in data.matrix[1:]]

        # Both ligands show up as columns, annotated with their pIC50 from mini_activity.csv
        assert any("7.8" in col for col in header)
        assert any("5.1" in col for col in header)

        # Both interacting residues show up as rows, sorted by residue number (41 < 142 < 166)
        assert row_labels == ["HIS 41", "ASN 142", "GLU 166"]

    def test_happy_path_with_template_file(self, analyzer, arpeggio_dir, activity_csv, template_path):
        data = analyzer.analyze_files(
            directory=arpeggio_dir,
            mode=analyzer.ARPEGGIO,
            activity_file=activity_csv,
            template_file=template_path,
        )
        # Template restricts interaction labels to those actually used in the ARPEGGIO_CONT/TYPE ordering
        assert "hbond" in data.interactions
        assert "plane-plane" in data.interactions

    def test_missing_activity_file_raises_file_not_found(self, analyzer, arpeggio_dir):
        with pytest.raises(FileNotFoundError):
            analyzer.analyze_files(
                directory=arpeggio_dir,
                mode=analyzer.ARPEGGIO,
                activity_file="does_not_exist.csv",
            )

    def test_missing_directory_raises(self, analyzer):
        with pytest.raises(FileOrDirectoryException):
            analyzer.analyze_files(directory="no_such_directory", mode=analyzer.ARPEGGIO)

    def test_invalid_mode_raises(self, analyzer, arpeggio_dir):
        with pytest.raises(InvalidModeException):
            analyzer.analyze_files(directory=arpeggio_dir, mode="not_a_real_mode")

    def test_save_writes_xlsx(self, analyzer, arpeggio_dir, activity_csv, tmp_path):
        analyzer.analyze_files(
            directory=arpeggio_dir,
            mode=analyzer.ARPEGGIO,
            activity_file=activity_csv,
            save="baseline.xlsx",
        )
        assert os.path.isfile(os.path.join(str(tmp_path), "baseline.xlsx"))

    def test_no_activity_file_labels_ligands_without_activity(self, analyzer, arpeggio_dir):
        data = analyzer.analyze_files(directory=arpeggio_dir, mode=analyzer.ARPEGGIO)
        header = data.matrix[0]
        # No activity file provided: ligand columns are just the ligand code, no "(value)" suffix
        assert all("(" not in col for col in header[1:])
