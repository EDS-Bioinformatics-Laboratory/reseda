from nose import with_setup
from nose.tools import assert_raises
from TranslateAndExtractCdr3 import getImgtMotifs
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
motif = getImgtMotifs()
p = regex.compile(motif, regex.BESTMATCH)


def test_no_match_results_in_none():
    ''' TranslateAndExtractCdr3.extractCDR3: No CDR3 found should result in cdr3pep None. Case: one mismatch too many'''
    peptide = "xxxxxxxxxxxxxATHabdeghVTVSyyyyyyyyy"
    (cdr3pep, aa_pos) = extractCDR3(cellType, peptide, p)
    assert cdr3pep == None
    assert aa_pos == []


def test_succesful_cdr3_extraction_cys_val():
    ''' TranslateAndExtractCdr3.extractCDR3: CDR3 recognizes a Cys and Val '''
    peptide = "xxxxxxxxxxxxxDTATHabCdeghVTVSyyyyyyyyy"
    (cdr3pep, aa_pos) = extractCDR3(cellType, peptide, p)
    assert cdr3pep == "CdeghV"
    assert aa_pos == [20, 26]


def test_succesful_cdr3_extraction_cys_trp():
    ''' TranslateAndExtractCdr3.extractCDR3: CDR3 recognizes a Cys and Trp '''
    peptide = "xxxxxxxxxxxxxDTATHabCdeghWyyyyyyyyy"
    (cdr3pep, aa_pos) = extractCDR3(cellType, peptide, p)
    assert cdr3pep == "CdeghW"
    assert aa_pos == [20, 26]


def test_succesful_cdr3_extraction_cys_phe():
    ''' TranslateAndExtractCdr3.extractCDR3: CDR3 recognizes a Cys and Phe '''
    peptide = "xxxxxxxxxxxxxDTATHabCdeghFyyyyyyyyy"
    (cdr3pep, aa_pos) = extractCDR3(cellType, peptide, p)
    assert cdr3pep == "CdeghF"
    assert aa_pos == [20, 26]


def test_multiple_motifs_shortest_match_is_preferred1():
    ''' TranslateAndExtractCdr3.extractCDR3: When multiple motifs are found the shortest should be returned'''
    peptide = "xxxxxxxCxxxxxDTATHabCdeghFWyyyyyyyyy"
    (cdr3pep, aa_pos) = extractCDR3(cellType, peptide, p)
    assert cdr3pep == "CdeghF"
    assert aa_pos == [20, 26]
