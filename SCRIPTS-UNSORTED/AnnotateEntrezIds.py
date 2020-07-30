from __future__ import print_function
import sys
import json

from Bio import Entrez

# *Always* tell NCBI who you are
Entrez.email = "b.d.vanschaik@amsterdamumc.nl"

def retrieve_annotation(id_list):

    """Annotates Entrez Gene IDs using Bio.Entrez, in particular epost (to
    submit the data to NCBI) and esummary to retrieve the information.
    Returns a list of dictionaries with the annotations."""

    request = Entrez.epost("gene",id=",".join(id_list))
    try:
        result = Entrez.read(request)
    except RuntimeError as e:
        #FIXME: How generate NAs instead of causing an error with invalid IDs?
        print("An error occurred while retrieving the annotations.")
        print("The error returned was %s" % e)
        sys.exit(-1)

    webEnv = result["WebEnv"]
    queryKey = result["QueryKey"]
    data = Entrez.esummary(db="gene", webenv=webEnv, query_key =
            queryKey)
    annotations = Entrez.read(data)

    print("Retrieved %d annotations for %d genes" % (len(annotations),
            len(id_list)))

    return(annotations)

myfile = "S233-head1000.tabular"
myids = list()
fh = open(myfile)
for line in fh:
    line = line.rstrip()
    myids.append(line)
fh.close()

fhOut = open(myfile + ".annotated", "w")
for myid in myids:
    js = retrieve_annotation(myid)
    print(myid, js["DocumentSummarySet"]["DocumentSummary"][0]["Name"], js["DocumentSummarySet"]["DocumentSummary"][0]["Description"], sep="\t", file=fhOut)
fhOut.close()
