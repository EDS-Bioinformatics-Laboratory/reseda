'''
@helpdesk: SURFsara helpdesk <helpdesk@surfsara.nl>

usage: python pilot.py
description:
    Connect to PiCaS server
    Get the next token in todo View
    Fetch the token parameters, e.g. input value
    Run main job (execute-all.sh) with the input argument
    When done, return the exit code to the token
    Attach the logs to the token
'''

#python imports
from __future__ import print_function
import os
import sys
import time
import shutil
import couchdb
import picasconfig
import json

#picas imports
from picas.actors import RunActor
from picas.clients import CouchClient
from picas.iterators import BasicViewIterator
from picas.modifiers import BasicTokenModifier
from picas.executers import execute

class ExampleActor(RunActor):
    def __init__(self, iterator, modifier):
        self.iterator = iterator
        self.modifier = modifier
        self.client = iterator.client

    def process_token(self, key, token):
        # Print token information
        print("-----------------------")
        print("Working on token: " + token['_id'])
        for k, v in token.iteritems():
            print(k, v)
            print("-----------------------")

        # Get the parameters
        js = json.loads(token['input'])

        # Write SAMPLES file
        try:
            fhOut = open("SAMPLES", "w")
        except:
            sys.exit("cannot write SAMPLES file")
        for sample in js["samples"]:
            print(sample, file=fhOut)
        fhOut.close()

        # Start running the main job
        command = "/usr/bin/time -v " + " ".join(["./execute-all.sh", "-r", js["run"], "-m", js["mids"], "-org", js["organism"], "-cell", js["cell"], "-celltype", js["celltype"], "-mm", str(js["mismatches"]), "-cregion", js["cregion"], "-p", js["protocol"], "-o", js["outdir"], "-b", js["barcodes"], "-u", js["umis"]])

        print(command)

        out = execute(command,shell=True)

        # Get the job exit code in the token
        token['exit_code'] = out[0]

        token = self.modifier.close(token)
        self.client.db[token['_id']] = token

        # Attach logs in token
        curdate=time.strftime("%d/%m/%Y_%H:%M:%S_")
        try:
            logsout = "logs_" + str(token['_id']) + ".out"
            log_handle = open(logsout, 'rb')
            self.client.db.put_attachment(token,log_handle,curdate+logsout)

            logserr = "logs_" + str(token['_id']) + ".err"
            log_handle = open(logserr, 'rb')
            self.client.db.put_attachment(token,log_handle,curdate+logserr)

        except:
            pass


    def cleanup_run(self):

        '''
        Description: Clean up input and result files
        In: -
        Out: -
        '''

        myfiles = os.listdir(".")
        for myfile in myfiles:
            if "L001" in myfile:
                if os.path.isfile(myfile):
                    os.remove(myfile)
                else:
                    shutil.rmtree(myfile, ignore_errors=True)

        myfiles = os.listdir("./split")
        for myfile in myfiles:
            if "L001" in myfile:
                shutil.rmtree("split/" + myfile, ignore_errors=True)

        myfiles = os.listdir("./split")
        for myfile in myfiles:
            if "L001" in myfile:
                os.remove("split/" + myfile)

        myfiles = os.listdir("./final")
        for myfile in myfiles:
            if "L001" in myfile or "mutations" in myfile:
                os.remove("final/" + myfile)


def main():
    # setup connection to db
    client = CouchClient(url=picasconfig.PICAS_HOST_URL, db=picasconfig.PICAS_DATABASE, username=picasconfig.PICAS_USERNAME, password=picasconfig.PICAS_PASSWORD)
    # Create token modifier
    modifier = BasicTokenModifier()
    # Create iterator, point to the right todo view
    iterator = BasicViewIterator(client, "Monitor/todo", modifier)
    # Create actor
    actor = ExampleActor(iterator, modifier)
    # Start work!
    print("Connected to the database %s sucessfully. Now starting work..." %(picasconfig.PICAS_DATABASE))
    actor.run()

if __name__ == '__main__':
    main()
