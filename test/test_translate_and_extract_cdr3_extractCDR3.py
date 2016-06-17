from nose import with_setup
from nose.tools import assert_raises
from TranslateAndExtractCdr3 import extractCDR3
import regex

'''Tests for extractCDR3 in TranslateAndExtractCdr3.py'''

cellType = "IGH_HUMAN"
motif_exact = ".+(DTATH.+?VTVS)"
motif_inexact = ".+(DTATH.+?VTVS){e<=1}"
p_exact = regex.compile(motif_exact, regex.BESTMATCH)
p_inexact = regex.compile(motif_inexact, regex.BESTMATCH)

def test_succesful_cdr3_extraction_exact_match():
    ''' CDR3 recognized with exact match and a Cys '''
    peptide = "xxxxxxxxxxxxxDTATHabCdefgVTVSyyyyyyyyy"
    (cdr3pep, aa_pos) = extractCDR3(cellType, peptide, p_exact)
    assert cdr3pep == "CdefgVT"
    assert aa_pos == [20,27]

def test_succesful_cdr3_extraction_exact_match_no_cys():
    ''' CDR3 recognized with exact match and without a Cys '''
    peptide = "xxxxxxxxxxxxxDTATHxabcdefgVTVSyyyyyyyyy"
    (cdr3pep, aa_pos) = extractCDR3(cellType, peptide, p_exact)
    assert cdr3pep == "abcdefgVT"
    assert aa_pos == [19,28]

def test_succesful_cdr3_extraction_inexact_match():
    ''' CDR3 recognized with an inexact match and a Cys '''
    peptide = "xxxxxxxxxxxxxDzATHabCdefgVTVSyyyyyyyyy"
    (cdr3pep, aa_pos) = extractCDR3(cellType, peptide, p_inexact)
    assert cdr3pep == "CdefgVT"
    assert aa_pos == [20,27]

def test_succesful_cdr3_extraction_inexact_match_no_cys():
    ''' CDR3 recognized with an inexact match and without a Cys '''
    peptide = "xxxxxxxxxxxxxDzATHxabcdefgVTVSyyyyyyyyy"
    (cdr3pep, aa_pos) = extractCDR3(cellType, peptide, p_inexact)
    assert cdr3pep == "abcdefgVT"
    assert aa_pos == [19,28]

def test_succesful_cdr3_extraction_inexact_match_with_exact_patter():
    ''' In-exact match should also recognize CDR3 with exact match'''
    peptide = "xxxxxxxxxxxxxDTATHabCdefgVTVSyyyyyyyyy"
    (cdr3pep, aa_pos) = extractCDR3(cellType, peptide, p_inexact)
    assert cdr3pep == "CdefgVT"
    assert aa_pos == [20,27]

def test_no_match_results_in_none():
    ''' No CDR3 found should result in cdr3pep None. Case: one mismatch too many'''
    peptide = "xxxxxxxxxxxxxATHabCdefgVTVSyyyyyyyyy"
    (cdr3pep, aa_pos) = extractCDR3(cellType, peptide, p_inexact)
    assert cdr3pep == None
    assert aa_pos == []

def test_multiple_motifs_found_exact_without_cys():
    ''' When multiple motifs are found, the second should be the start of the CDR3 '''
    peptide = "xxxxxxxxxxxxxDTATHDTATHxabcdefgVTVSyyyyyyyyy"
    (cdr3pep, aa_pos) = extractCDR3(cellType, peptide, p_exact)
    assert cdr3pep == "abcdefgVT"
    assert aa_pos == [24,33]

def test_multiple_motifs_inexact_match_exact_is_preferred1():
    ''' When multiple motifs are found with inexact matching, exact match is preferred '''
    peptide = "xxxxxxxxxxxxxDTATHDzATHxabcdefgVTVSyyyyyyyyy"
    (cdr3pep, aa_pos) = extractCDR3(cellType, peptide, p_inexact)
    assert cdr3pep == "zATHxabcdefgVT"
    assert aa_pos == [19,33]

def test_multiple_motifs_inexact_match_exact_is_preferred2():
    ''' When multiple motifs are found with inexact matching, exact match is preferred '''
    peptide = "xxxxxxxxxxxxxDzATHDTATHxabcdefgVTVSyyyyyyyyy"
    (cdr3pep, aa_pos) = extractCDR3(cellType, peptide, p_inexact)
    assert cdr3pep == "abcdefgVT"
    assert aa_pos == [24,33]
