from nose.tools import assert_raises
from CountFirstTenNucleotides import getHighest

'''Tests for getHighest method in CountFirstTenNucleotides.py'''

def test_empty_dictionary():
    '''CountFirstTenNucleotides.getHighest: empty dictionary results in IndexError'''
    d = dict()
    assert_raises(IndexError, getHighest, d)

def test_good_dictionary():
    '''CountFirstTenNucleotides.getHighest: normal good dictionary returns highest key, value pair'''
    d = {"a":1000, "b":500, "c":20000, "d":2}
    assert getHighest(d) == ("c",20000)

def test_multiple_highest_values():
    '''CountFirstTenNucleotides.getHighest: dictionary with multiple highest values returns one of them at random'''
    d = {"a":20000, "b":500, "c":20000, "d":2}
    assert getHighest(d) == ("c",20000) or getHighest(d) == ("a", 20000)
