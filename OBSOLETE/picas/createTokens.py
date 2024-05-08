'''
@helpdesk: SURFsara helpdesk <helpdesk@surfsara.nl>

usage: python createTokens.py [path to tokens file]
description:
   Connects to PiCaS server
   Creates one token for each line in [tokens file]
   Loads the tokens
'''

import sys
import os
import couchdb
import random
import picasconfig

def loadTokens(db):
    tokens = []
    tokensfiles = sys.argv[1:]

    i = 0
    for tokenfile in tokensfiles:
        fh = open(tokenfile)
        input = fh.read().rstrip()
        token = {
            '_id': 'token_' + str(i),
            'type': 'token',
            'lock': 0,
            'done': 0,
            'hostname': '',
            'scrub_count': 0,
            'input': input,
            'exit_code': ''
        }
        tokens.append(token)
        i = i +1
        fh.close()
    db.update(tokens)

def get_db():
    server = couchdb.Server(picasconfig.PICAS_HOST_URL)
    username = picasconfig.PICAS_USERNAME
    pwd = picasconfig.PICAS_PASSWORD
    server.resource.credentials = (username,pwd)
    db = server[picasconfig.PICAS_DATABASE]
    return db

if __name__ == '__main__':
   if len(sys.argv) < 2:
       sys.exit("Usage: python createTokens.py ../tokens/*")
   #Create a connection to the server
   db = get_db()
   #Load the tokens to the database
   loadTokens(db)
