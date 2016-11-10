from __future__ import print_function
import sys
import csv
import json
import pymongo

def parseMiseqSheet (f):
    '''
    Description: Read sample sheet and put it in json format
    In: filename
    Out: json object
    '''

    js = dict()
    section = "doc"
    keys = list()

    with open(f) as csvfile:
        fh = csv.reader(csvfile)
        for c in fh:

            # print(c)

            # Which section?
            if len(c) == 0:   # skip empty
                continue
            elif c[0] == "":   # skip empty
                continue
            elif c[0] == "[Header]":
                section = "Header"
                continue
            elif c[0] == "[Reads]":
                section = "Reads"
                continue
            elif c[0] == "[Settings]":
                section = "Settings"
                continue
            elif c[0] == "[Data]":
                section = "Data"
                continue

            if section == "Header" or section == "Settings":
                js[c[0]] = c[1]  # first column is key, second is the value
            elif section == "Data":
                if len(keys) == 0:  # column names -> keys
                    keys = c[:]
                    js["Samples"] = list()
                else:               # data, columns are values
                    d = dict()
                    for i in range(len(c)):  # first column (c[0]) is the sample name, rest is meta data
                        d[keys[i]] = c[i]
                    js["Samples"].append(d)

    csvfile.close()
    return(js)

if __name__ == '__main__':

    if len(sys.argv) < 2:
        sys.exit("Usage: " + sys.argv[0] + " Miseq-sample-sheet.csv")

    samplesheet = sys.argv[1]

    # client = pymongo.MongoClient()
    # db = client.immuno

    js = parseMiseqSheet(samplesheet)
    print(json.dumps(js, indent=4))

    # result = db.samples.insert_one(js)
    # print("Data stored:", result.inserted_id)
