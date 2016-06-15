import unittest
from sequences import complement

class SequencesTestCase(unittest.TestCase):
    """Tests for `sequences.py`."""

    def test_complement_convert_c_to_g(self):
        """Is CGAT successfully converted to GCTA?"""
        self.assertEqual(complement("CGAT"), "GCTA")

    def test_complement_convert_lowercase_letters(self):
    	''' does the method take complement of lowercase letters (should be "no")'''
    	self.assertEqual(complement("cgat"), "gcta")

if __name__ == '__main__':
    unittest.main()
