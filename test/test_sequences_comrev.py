from nose.tools import assert_raises
from sequences import comrev

'''Tests for complement method in sequences.py'''


def test_comrev_reverse_string():
    '''sequences.comrev: Is CGATM successfully converted to KATCG?'''
    assert comrev("CGATM") == "KATCG"


def test_comrev_empty_string():
    '''sequences.comrev: Empty string should return an empty string'''
    assert comrev("") == ""


def test_comrev_unknown_letters():
    '''sequences.comrev: Unknown letters should raise KeyError'''
    assert_raises(KeyError, comrev, "CATGXYZ")
