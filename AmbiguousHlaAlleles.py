from __future__ import print_function
import xml.etree.ElementTree as ETree

def getLookupTable (f):
    '''
    Description: Read XML file with ambiguous HLA names and store in a dictionary
    In: f (filename)
    Out: lookupTable[ambig_allele] = main_allele (dict)
    '''
   
    print("=== Reading XML file ===")
    tree = ETree.ElementTree(file=hla_ambig_file)
    print("=== done ===")

    lookupTable = dict()

    # Get main group and the ambiguous allele names
    for node in tree.findall('.//{http://www.example.org/ambig-aw}gGroup'):
        main_allele = node.attrib['name']
        # Get child nodes
        for child in node.getchildren():
            ambig_allele = child.attrib['name']
            lookupTable[ambig_allele] = main_allele
            print(ambig_allele, main_allele)

    return(lookupTable)

if __name__ == '__main__':
    hla_ambig_file = "/home/barbera/git/IMGTHLA/xml/test_hla_ambigs.xml"
    lookupTable = getLookupTable(hla_ambig_file)
