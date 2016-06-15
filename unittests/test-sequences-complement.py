from nose.tools import assert_raises
from sequences import complement

'''Tests for complement method in sequences.py'''

def test_complement_convert_c_to_g():
    '''Is CGATM successfully converted to GCTAK?'''
    assert complement("CGATM") == "GCTAK"

def test_complement_empty_string():
    '''Complement of empty string should return an empty string'''
    assert complement("") == ""

def test_complement_keyerror_lowercase_sequence():
    '''Method doesn't take complement of lowercase letters'''
    assert_raises(KeyError, complement, "cgat")

def test_complement_keyerror_general():
    ''' Raise a KeyError when an unknown letter is provided'''
    assert_raises(KeyError, complement, "XYZCATG")

def test_complement_keyerror_something_else_then_string():
    ''' Raise a TypeError when a number is provided (not an iterable object)'''
    assert_raises(TypeError, complement, 123)
