import unittest
from sequences import complement

class SequencesTestCase(unittest.TestCase):
    '''Tests for complement method in sequences.py'''

    def test_complement_convert_c_to_g(self):
        '''Is CGATM successfully converted to GCTAK?'''
        self.assertEqual(complement("CGATM"), "GCTAK")

    def test_complement_keyerror_lowercase_sequence(self):
    	''' does the method take complement of lowercase letters (should be "no")'''
        with self.assertRaises(KeyError):
    	   complement("cgat")

    def test_complement_empty_string(self):
        '''Should return an empty string'''
        self.assertEqual(complement(""), "")

    def test_complement_keyerror_general(self):
        ''' does the method raise a keyerror when an unknown letter is provided'''
        with self.assertRaises(KeyError):
           complement("XYZCATG")

    def test_complement_keyerror_something_else_then_string(self):
        ''' does the method raise a TypeError when a number is provided'''
        with self.assertRaises(TypeError):
           complement(123)

if __name__ == '__main__':
    unittest.main()
