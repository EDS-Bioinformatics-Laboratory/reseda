from nose import with_setup
from nose.tools import assert_raises
from TranslateAndExtractCdr3 import extractCDR3
import regex

''' Tests for extractCDR3 in TranslateAndExtractCdr3.py '''

'''
CDR3 found / not found
CDR3's with and without Cys
Exact matching and in-exact matching
In-exact matching should find exact match
Multiple matches: two exact matches (take most downstream one), exact match+inexact match and in-exact match+exact match (take best match)
Different cell types (not tested yet)
Search with multiple motifs (not tested yet)
'''

cellType = "IGH_HUMAN"
motif_exact = ".+(DTATH.+?VTVS)"
motif_inexact = ".+(DTATH.+?VTVS){e<=1}"
motif_full_inexact = ".+(DTATH.+?VTVS|GTAVY.+?VTVS|DMTVY.+?VTVS|GTVVY.+?VTVS|GTAAY.+?VTVS|DMAVY.+?VTVS|DVAVY.+?VTVS|DAAMY.+?VTVS|DTATY.+?VTVS|DMTMH.+?VTVS|DTVVY.+?VTVS|DTALY.+?VTVS|DSAVY.+?VTVS|DTAVY.+?VTVS|DMAMY.+?VTVS|DTAMY.+?VTVS){e<=1}"
p_exact = regex.compile(motif_exact, regex.BESTMATCH)
p_inexact = regex.compile(motif_inexact, regex.BESTMATCH)
p_full_inexact = regex.compile(motif_full_inexact, regex.BESTMATCH)

def test_no_match_results_in_none():
    ''' TranslateAndExtractCdr3.extractCDR3: No CDR3 found should result in cdr3pep None. Case: one mismatch too many'''
    peptide = "xxxxxxxxxxxxxATHabCdefgVTVSyyyyyyyyy"
    (cdr3pep, aa_pos) = extractCDR3(cellType, peptide, p_inexact)
    assert cdr3pep == None
    assert aa_pos == []

def test_succesful_cdr3_extraction_exact_match():
    ''' TranslateAndExtractCdr3.extractCDR3: CDR3 recognized with exact match and a Cys '''
    peptide = "xxxxxxxxxxxxxDTATHabCdefgVTVSyyyyyyyyy"
    (cdr3pep, aa_pos) = extractCDR3(cellType, peptide, p_exact)
    assert cdr3pep == "CdefgVT"
    assert aa_pos == [20,27]

def test_succesful_cdr3_extraction_exact_match_no_cys():
    ''' TranslateAndExtractCdr3.extractCDR3: CDR3 recognized with exact match and without a Cys '''
    peptide = "xxxxxxxxxxxxxDTATHxabcdefgVTVSyyyyyyyyy"
    (cdr3pep, aa_pos) = extractCDR3(cellType, peptide, p_exact)
    assert cdr3pep == "abcdefgVT"
    assert aa_pos == [19,28]

def test_succesful_cdr3_extraction_inexact_match():
    ''' TranslateAndExtractCdr3.extractCDR3: CDR3 recognized with an inexact match and a Cys '''
    peptide = "xxxxxxxxxxxxxDzATHabCdefgVTVSyyyyyyyyy"
    (cdr3pep, aa_pos) = extractCDR3(cellType, peptide, p_inexact)
    assert cdr3pep == "CdefgVT"
    assert aa_pos == [20,27]

def test_succesful_cdr3_extraction_inexact_match_no_cys():
    ''' TranslateAndExtractCdr3.extractCDR3: CDR3 recognized with an inexact match and without a Cys '''
    peptide = "xxxxxxxxxxxxxDzATHxabcdefgVTVSyyyyyyyyy"
    (cdr3pep, aa_pos) = extractCDR3(cellType, peptide, p_inexact)
    assert cdr3pep == "abcdefgVT"
    assert aa_pos == [19,28]

def test_succesful_cdr3_extraction_inexact_match_with_exact_patter():
    ''' TranslateAndExtractCdr3.extractCDR3: In-exact match should also recognize CDR3 with exact match'''
    peptide = "xxxxxxxxxxxxxDTATHabCdefgVTVSyyyyyyyyy"
    (cdr3pep, aa_pos) = extractCDR3(cellType, peptide, p_inexact)
    assert cdr3pep == "CdefgVT"
    assert aa_pos == [20,27]

def test_multiple_motifs_found_exact_without_cys():
    ''' TranslateAndExtractCdr3.extractCDR3: When multiple motifs are found, the last one should be the start of the CDR3 '''
    peptide = "xxxxxxxxxxxxxDTATHDTATHxabcdefgVTVSyyyyyyyyy"
    (cdr3pep, aa_pos) = extractCDR3(cellType, peptide, p_exact)
    assert cdr3pep == "abcdefgVT"
    assert aa_pos == [24,33]

def test_multiple_motifs_inexact_match_exact_is_preferred1():
    ''' TranslateAndExtractCdr3.extractCDR3: When multiple motifs are found with inexact matching, exact match is preferred (1) '''
    peptide = "xxxxxxxxxxxxxDTATHDzATHxabcdefgVTVSyyyyyyyyy"
    (cdr3pep, aa_pos) = extractCDR3(cellType, peptide, p_inexact)
    assert cdr3pep == "zATHxabcdefgVT"
    assert aa_pos == [19,33]

def test_multiple_motifs_inexact_match_exact_is_preferred2():
    ''' TranslateAndExtractCdr3.extractCDR3: When multiple motifs are found with inexact matching, exact match is preferred (2) '''
    peptide = "xxxxxxxxxxxxxDzATHDTATHxabcdefgVTVSyyyyyyyyy"
    (cdr3pep, aa_pos) = extractCDR3(cellType, peptide, p_inexact)
    assert cdr3pep == "abcdefgVT"
    assert aa_pos == [24,33]
