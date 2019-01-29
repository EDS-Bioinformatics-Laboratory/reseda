from nose.tools import assert_raises
from TranslateAndExtractCdr3 import getMotifs

'''Tests for getMotifs in TranslateAndExtractCdr3.py'''


def test_getMotifs_list_exact_match():
    '''TranslateAndExtractCdr3.getMotifs: Do we get the motifs with pattern for exact match?'''
    assert getMotifs("IGH_HUMAN", 0).endswith('VTVS)')

def test_getMotifs_list_1_match():
    '''TranslateAndExtractCdr3.getMotifs: Do we get the motifs with 1 mismatch?'''
    assert getMotifs("IGH_HUMAN", 1).endswith('VTVS){e<=1}')

def test_getMotifs_error_unknown_celltype():
    ''' TranslateAndExtractCdr3.getMotifs: Trying to get the motifs for an unknown cell type should result in an error'''
    assert_raises(NameError, getMotifs, "blah", 0)


def test_getMotifs_error_negative_mismatch():
    '''TranslateAndExtractCdr3.getMotifs: Defining a negative number for mismatch should result in an error'''
    assert_raises(TypeError, getMotifs, "IGH_HUMAN", -1)


def test_getMotifs_check_error_wrong_type():
    '''TranslateAndExtractCdr3.getMotifs: Defining any other type than 0 or a positive integer should result in an error'''
    assert_raises(TypeError, getMotifs, "IGH_HUMAN", "a")
