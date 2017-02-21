import openpyxl as xlsx
import xlrd as xls
import sys


def findInRow(p, s):
    '''
    In: p = search terms (list), s = string
    Out: True or False (found a match or not)
    '''
    match = False
    for pattern in p:
        if pattern in s:
            match = True

    return(match)


def readXlsx(f, p):
    # DOC: https://openpyxl.readthedocs.io/en/default/tutorial.html#loading-from-a-file
    wb = xlsx.load_workbook(f)
    sheets = wb.get_sheet_names()
    print(sheets)
    for s in sheets:
        sheet_ranges = wb[s]
        print(s, sheet_ranges.max_row, 'x', sheet_ranges.max_column)


def readXls(f, p):
    # DOC: https://github.com/python-excel/tutorial/blob/master/python-excel.pdf
    wb = xls.open_workbook(filename=f)
    sheets = list()
    for s in wb.sheets():
        sheets.append(s.name)
        print(s.name, s.nrows, 'x', s.ncols)

        for row in range(s.nrows):
            values = list()
            for col in range(s.ncols):
                values.append(str(s.cell(row, col).value))

        myrow = ','.join(values)
        if findInRow(p, myrow) is True:
            print(myrow)

    print(sheets)


if __name__ == '__main__':
    # Check extension: .xls .xlsx .csv and use the appropriate package to open the file
    f = sys.argv[1]
    ext = f.split(".")[-1]

    p = list()
    for i in range(25):
        p.append("HD" + str(i))
    print("p = ", p)

    if ext == "xlsx":
        readXlsx(f, p)
    elif ext == "xls":
        readXls(f, p)
    else:
        print("Unknown extension:", ext)

    # Open each sheet and search for a particular term (or terms)
    # Print the entire row if search term is found
