from nose.tools import assert_raises
from nose.tools import assert_dict_equal
from CountFirstTenNucleotides import countFirstNucleotides

'''Tests for countFirstNucleotides method in CountFirstTenNucleotides.py'''

def test_good_example():
    '''CountFirstTenNucleotides.countFirstNucleotides: good example should go well'''
    myfile = "test/data/testdata_CountFirstTenNucleotides_tabdel_second_column.csv"
    assert_dict_equal(countFirstNucleotides(myfile,2,1), {"aa":3,"bb":2,"cc":1,"dd":1,"ee":1})

def test_data_in_last_column():
    '''CountFirstTenNucleotides.countFirstNucleotides: get data from last column'''
    myfile = "test/data/testdata_CountFirstTenNucleotides_last_column.csv"
    assert_dict_equal(countFirstNucleotides(myfile,2,1), {"aa":3,"bb":2,"cc":1,"dd":1,"ee":1})

def test_not_tab_delimited_gives_warning():
    '''CountFirstTenNucleotides.countFirstNucleotides: not tab delimited, but space delimited, raises warning'''
    myfile = "test/data/testdata_CountFirstTenNucleotides_not_tab_delimited.csv"
    assert_raises(Warning, countFirstNucleotides, myfile, 2, 1)

def test_not_tab_delimited_get_first_col_warning():
    '''CountFirstTenNucleotides.countFirstNucleotides: not tab delimited, but space delimited and getting first column gives a warning'''
    myfile = "test/data/testdata_CountFirstTenNucleotides_not_tab_delimited.csv"
    assert_raises(Warning, countFirstNucleotides, myfile, 2, 0)

def test_file_with_one_column():
    '''CountFirstTenNucleotides.countFirstNucleotides: file with just one column should go well'''
    myfile = "test/data/testdata_CountFirstTenNucleotides_one_column.csv"
    assert_dict_equal(countFirstNucleotides(myfile,2,0), {"aa":3,"bb":2,"cc":1,"dd":1,"ee":1})
