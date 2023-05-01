# Code from: https://www.biostars.org/p/15011/

import sys, operator

header = ''
data = ''

allData = []

for line in sys.stdin:
    if line[0] == "@":
        if data != '':
            allData.append((header,data))

        header = line.strip()[1:]
        data = line
    else:
        data += line

allData.append((header,data))
allData.sort(key = operator.itemgetter(0))

for item in allData:
    print(item[1], end='')
