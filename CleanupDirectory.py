import os
import shutil


def cleanUp():
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

    if os.path.isdir("./orig"):
        myfiles = os.listdir("./orig")
        for myfile in myfiles:
            if os.path.isfile("orig/" + myfile):
                os.remove("orig/" + myfile)
            else:
                shutil.rmtree("orig/" + myfile, ignore_errors=True)

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
            if os.path.isfile("final/" + myfile):
                os.remove("final/" + myfile)
            else:
                shutil.rmtree("final/" + myfile, ignore_errors=True)


if __name__ == '__main__':
    cleanUp()
