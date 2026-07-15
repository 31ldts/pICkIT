"""Tests for pickit.parsers (common helpers, arpeggio, ichem).

New tests introduced during Fase 3, when the parsers were split out of
io_mixin.py. Arpeggio behavior is already exercised end-to-end through
test_analyze_files.py (Fase 0 baseline, frozen); these tests target the
parser functions directly and in isolation, including an IChem case using
hand-built line content — the Fase 0 notes explicitly flagged that no
IChem sample data was available among the provided files, so
`parse_ichem_file` had 0% dedicated coverage until now.
"""

from pickit.parsers.common import get_protein_ligand, modify_cell, validate_string
from pickit.parsers.ichem import parse_ichem_file


class TestModifyCell:
    def test_empty_cell_adds_first_interaction(self):
        result = modify_cell(text="", interaction="hbond", atoms="CA-O1(A)", interaction_labels=["hbond", "ionic"])
        assert result == "1 |CA-O1(A)|"

    def test_appends_new_interaction_type_to_existing_cell(self):
        result = modify_cell(
            text="1 |CA-O1(A)|", interaction="ionic", atoms="CB-N1(A)", interaction_labels=["hbond", "ionic"]
        )
        assert result == "1 |CA-O1(A)|; 2 |CB-N1(A)|"

    def test_adds_atoms_to_existing_interaction_type(self):
        result = modify_cell(
            text="1 |CA-O1(A)|", interaction="hbond", atoms="CB-N1(A)", interaction_labels=["hbond", "ionic"]
        )
        assert result == "1 |CA-O1(A), CB-N1(A)|"

    def test_does_not_duplicate_identical_atom_pair(self):
        result = modify_cell(
            text="1 |CA-O1(A)|", interaction="hbond", atoms="CA-O1(A)", interaction_labels=["hbond", "ionic"]
        )
        assert result == "1 |CA-O1(A)|"

    def test_unknown_interaction_falls_back_to_last_label(self):
        result = modify_cell(
            text="", interaction="not_a_label", atoms="CA-O1(A)", interaction_labels=["hbond", "ionic"]
        )
        assert result == "2 |CA-O1(A)|"


class TestValidateString:
    def test_accepts_well_formed_residue_string(self):
        assert validate_string("HIS 41-A") is True

    def test_rejects_missing_dash(self):
        assert validate_string("HIS 41") is False

    def test_rejects_lowercase_residue_code(self):
        assert validate_string("his 41-A") is False


class TestGetProteinLigand:
    AA = ["HIS", "ASN", "GLU"]

    def test_returns_protein_then_ligand_when_only_begin_is_protein(self):
        begin = {"label_comp_type": "P", "label_comp_id": "HIS"}
        end = {"label_comp_type": "N", "label_comp_id": "LIG"}
        prot, lig = get_protein_ligand(begin=begin, end=end, amino_acid_codes=self.AA)
        assert prot is begin and lig is end

    def test_returns_protein_then_ligand_when_only_end_is_protein(self):
        begin = {"label_comp_type": "N", "label_comp_id": "LIG"}
        end = {"label_comp_type": "P", "label_comp_id": "HIS"}
        prot, lig = get_protein_ligand(begin=begin, end=end, amino_acid_codes=self.AA)
        assert prot is end and lig is begin

    def test_both_polymer_disambiguates_by_amino_acid_code(self):
        begin = {"label_comp_type": "P", "label_comp_id": "XYZ"}  # not a standard amino acid
        end = {"label_comp_type": "P", "label_comp_id": "HIS"}
        prot, lig = get_protein_ligand(begin=begin, end=end, amino_acid_codes=self.AA)
        assert prot is end and lig is begin

    def test_both_polymer_and_both_standard_amino_acids_returns_none(self):
        begin = {"label_comp_type": "P", "label_comp_id": "HIS"}
        end = {"label_comp_type": "P", "label_comp_id": "ASN"}
        prot, lig = get_protein_ligand(begin=begin, end=end, amino_acid_codes=self.AA)
        assert prot is None and lig is None

    def test_neither_is_protein_returns_none(self):
        begin = {"label_comp_type": "N", "label_comp_id": "LIG"}
        end = {"label_comp_type": "N", "label_comp_id": "LG2"}
        prot, lig = get_protein_ligand(begin=begin, end=end, amino_acid_codes=self.AA)
        assert prot is None and lig is None


class TestParseIchemFile:
    def test_single_interaction_line_populates_matrix(self):
        # IChem line format: 10 '|'-separated fields; [0]=interaction,
        # [1]=protein atom, [3]=residue ("RES NUM-SUBUNIT"), [4]=ligand atom.
        # "Hydrophobic" is INTERACTION_LABELS[0] (IChem's own label set, not Arpeggio's).
        line = "Hydrophobic|CA|x|HIS 41-A|O1|x|x|x|x|x"
        matrix, aa, cont, subunits_set = parse_ichem_file(
            content=[line],
            index=0,
            files=["complex1.txt"],
            subunits_set=set(),
            cont=0,
            matrix=[],
            aa={},
            protein=True,
            ligand=True,
            subunit=False,
        )
        assert aa == {"HIS 41": 0}
        assert cont == 1
        assert subunits_set == {"A"}
        assert matrix[0][0] == "1 |CA-O1(A)|"

    def test_malformed_line_is_ignored(self):
        matrix, aa, cont, subunits_set = parse_ichem_file(
            content=["not|enough|fields"],
            index=0,
            files=["complex1.txt"],
            subunits_set=set(),
            cont=0,
            matrix=[],
            aa={},
            protein=True,
            ligand=True,
            subunit=False,
        )
        assert matrix == []
        assert aa == {}

    def test_subunit_true_keeps_subunit_in_residue_label(self):
        line = "Hydrophobic|CA|x|HIS 41-A|O1|x|x|x|x|x"
        matrix, aa, cont, subunits_set = parse_ichem_file(
            content=[line],
            index=0,
            files=["complex1.txt"],
            subunits_set=set(),
            cont=0,
            matrix=[],
            aa={},
            protein=True,
            ligand=True,
            subunit=True,
        )
        assert aa == {"HIS 41-A": 0}
        # subunits_set is only populated in the not-subunit branch
        assert subunits_set == set()
