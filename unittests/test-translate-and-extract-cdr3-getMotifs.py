from nose.tools import assert_raises
from TranslateAndExtractCdr3 import getMotifs

'''Tests for getMotifs in TranslateAndExtractCdr3.py'''

def test_getMotifs_list_exact_match():
    '''Do we get the motifs with pattern for exact match?'''
    assert getMotifs("IGH_HUMAN", 0) == '.+(DTATH.+?VTVS|GTAVY.+?VTVS|DMTVY.+?VTVS|GTVVY.+?VTVS|GTAAY.+?VTVS|DMAVY.+?VTVS|DVAVY.+?VTVS|DAAMY.+?VTVS|DTATY.+?VTVS|DMTMH.+?VTVS|DTVVY.+?VTVS|DTALY.+?VTVS|DSAVY.+?VTVS|DTAVY.+?VTVS|DMAMY.+?VTVS|DTAMY.+?VTVS)'

def test_getMotifs_list_1_match():
    '''Do we get the motifs with 1 mismatch?'''
    assert getMotifs("IGH_HUMAN", 1) == '.+(DTATH.+?VTVS|GTAVY.+?VTVS|DMTVY.+?VTVS|GTVVY.+?VTVS|GTAAY.+?VTVS|DMAVY.+?VTVS|DVAVY.+?VTVS|DAAMY.+?VTVS|DTATY.+?VTVS|DMTMH.+?VTVS|DTVVY.+?VTVS|DTALY.+?VTVS|DSAVY.+?VTVS|DTAVY.+?VTVS|DMAMY.+?VTVS|DTAMY.+?VTVS){e<=1}'

def test_getMotifs_error_unknown_celltype():
    ''' Trying to get the motifs for an unknown cell type should result in an error'''
    assert_raises(NameError, getMotifs, "blah", 0)

def test_getMotifs_error_negative_mismatch():
    '''Defining a negative number for mismatch should result in an error'''
    assert_raises(TypeError, getMotifs, "IGH_HUMAN", -1)

def test_getMotifs_check_error_string():
    '''Defining any other type than 0 or a positive integer should result in an error'''
    assert_raises(TypeError, getMotifs, "IGH_HUMAN", "a")
