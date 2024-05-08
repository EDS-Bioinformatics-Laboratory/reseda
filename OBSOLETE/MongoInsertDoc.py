from __future__ import print_function
import sys
import json
import pymongo


if __name__ == '__main__':
    if len(sys.argv) < 2:
        sys.exit("Usage: " + sys.argv[0] + " *.json")

    # Connect to mongodb
    client = pymongo.MongoClient()
    db = client.immuno

    myfiles = sys.argv[1:]

    # Read documents
    for myfile in myfiles:
        fh = open(myfile)
        text = fh.read()
        fh.close()
        js = json.loads(text)

        # Insert document
        result = db.samples.insert_one(js)
        print("Data stored:", result.inserted_id)
