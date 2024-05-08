import couchdb
import picasconfig

def get_db():
    server = couchdb.Server(picasconfig.PICAS_HOST_URL)
    username = picasconfig.PICAS_USERNAME
    pwd = picasconfig.PICAS_PASSWORD
    server.resource.credentials = (username,pwd)
    db = server[picasconfig.PICAS_DATABASE]
    return db

if __name__ == '__main__':
   #Create a connection to the server
   db = get_db()

   print("To do:", len([x for x in db.iterview("Monitor/todo", 1000000)]))
   print("Done:", len([x for x in db.iterview("Monitor/done", 1000000)]))
   print("Locked:", len([x for x in db.iterview("Monitor/locked", 1000000)]))
   print("Error:", len([x for x in db.iterview("Monitor/error", 1000000)]))
