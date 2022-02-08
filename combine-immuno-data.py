from __future__ import print_function
import sqlite3
import argparse
import sys


def create_table(name, colnames):
    '''
    Description: create a table in the database
    In: name (table name), colnames (list with column names)
    Out: query (string)
    '''
    query = "CREATE TABLE " + name + " ("
    query = query + ", ".join(colnames)
    query = query + ");"
    print(query)
    cur.execute(query)
    return(query)


def import_data(datafile, delim, table, colnames):
    '''
    Description: Reads a csv file and inserts the content in a table
    In: datafile (path), delim (e.g. ',' or ' '), table (name of the table in database)
    Out: query (string)
    '''
    fh = open(datafile, 'r')
    for row in fh:
        row = row.rstrip()
        row = row.split(delim)
        query = "INSERT INTO " + table + " ("
        query = query + ", ".join(colnames)
        query = query + ") VALUES ("
        query = query + ", ".join(len(colnames) * ["?"]) + ");"
        # print(query, row)
        cur.execute(query, row)
    con.commit()


def clean_name(gene):
    '''
    Description: cleans up unnessary info in the gene name that was returned by the aligner, also removes allele info
    In: V or J name (e.g. 'M99672|IGHV3-43*01|Homo')
    Out: cleaned up gene name, sub, and main (e.g. 'IGHV3-43*01', 'IGHV3-43', 'IGHV3')
    '''
    tmp = gene.split("|")

    if len(tmp) == 3:
        gene = tmp[1]   # remove accession at start and 'Homo' or 'Mouse' at the end
        gene_sub = gene.split("*")[0]   # remove allele info
        gene_main = gene_sub.split("-")[0]  # remove sub info
        return(gene, gene_sub, gene_main)
    else:
        return(gene, gene, gene)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Combine tables that were produced in earlier steps')
    parser.add_argument('-q', '--qual', default=30, type=int, help='Minimum base quality of CDR3 (default: %(default)s)')
    parser.add_argument('-m', '--midfile', default='SAMPLE_Sn_L001.assembled-report.txt', type=str, help='File with the identified MIDs (default: %(default)s)')
    parser.add_argument('-c', '--cdr3file', default='SAMPLE_Sn_L001.assembled.fastq.gz-IGH_HUMAN-CDR3.csv', type=str, help='File with identified CDR3s (default: %(default)s)')
    parser.add_argument('-v', '--vfile', default='SAMPLE_Sn_L001.assembled-IGHV_human-easy-import.txt', type=str, help='File with identified V regions (default: %(default)s)')
    parser.add_argument('-j', '--jfile', default='SAMPLE_Sn_L001.assembled-IGHJ_human-easy-import.txt', type=str, help='File with identified V regions (default: %(default)s)')
    parser.add_argument('-s', '--seqfile', default='SAMPLE_Sn_L001.assembled.fastq.gz-IGH_HUMAN.csv', type=str, help='File with ... to check ... (default: %(default)s)')
    parser.add_argument('-e', '--extrafile', default='SAMPLE_Sn_L001.assembled.fastq.gz-IGH_HUMAN-extra.txt', type=str, help='File with extra motifs, e.g. N-glycosylation sites (default: %(default)s)')
    parser.add_argument('-o', '--outfile', default='-all_info.csv', type=str, help='Name of the all_info output file (default: %(default)s)')
    parser.add_argument('-ocf', '--clonesfile', default='-clones.csv', type=str, help='Name of the clones output file (default: %(default)s)')
    parser.add_argument('-ocs', '--clonessubsfile', default='-clones-subs.csv', type=str, help='Name of the clones-subs output file (default: %(default)s)')
    parser.add_argument('-ocm', '--clonesmainsfile', default='-clones-mains.csv', type=str, help='Name of the clones-mains output file (default: %(default)s)')
    parser.add_argument('-t', '--totalfile', default='-productive.txt', type=str, help='Name of the summary (total) output file (default: %(default)s)')
    args = parser.parse_args()

    # Check if all arguments were specified
    for arg in vars(args):
        if str(getattr(args, arg)).startswith("SAMPLE_Sn") or str(getattr(args, arg)).startswith("-"):
            parser.print_help()
            exit()

    #################

    # Output file
    fhOut = open(args.outfile, 'w')
    fhClones = open(args.clonesfile, 'w')
    fhClonesSubs = open(args.clonessubsfile, 'w')
    fhClonesMains = open(args.clonesmainsfile, 'w')
    fhTotal = open(args.totalfile, 'w')


    con = sqlite3.connect(":memory:")
    # con = sqlite3.connect("test.db")
    cur = con.cursor()
    cur.execute("PRAGMA temp_store_directory='.'")

    # ********* Import all data ***********

    # MID
    colnames = ["acc", "beforeMID", "MID", "afterMID"]
    create_table("mid", colnames)
    import_data(args.midfile, " ", "mid", colnames)

    # CDR3
    colnames = ["acc", "readingframe", "cdr3pep", "cdr3nuc", "cdr3pepshort", "cdr3nucshort", "cdr3_qual_min", "cdr3_qual_max", "cdr3_qual_avg", "cdr3_qual", "nt_start", "nt_end", "seq_length"]
    create_table("cdr3", colnames)
    import_data(args.cdr3file, "\t", "cdr3", colnames)

    # V genes
    colnames = ["acc", "V_flag", "V_gene"]
    create_table("v", colnames)
    import_data(args.vfile, "\t", "v", colnames)

    # J genes
    colnames = ["acc", "J_flag", "J_gene"]
    create_table("j", colnames)
    import_data(args.jfile, "\t", "j", colnames)

    # Sequences
    colnames = ["acc", "readingframe_seq", "seq", "pep", "qual"]
    create_table("seq", colnames)
    import_data(args.seqfile, "\t", "seq", colnames)

    # N glycosylation sites
    colnames = ["acc", "readingframe", "site", "start", "end"]
    create_table("sites", colnames)
    import_data(args.extrafile, " ", "sites", colnames)

    # Combine all tables into one big table: all_info
    query = "CREATE TABLE all_info AS SELECT * FROM mid JOIN cdr3 USING (acc) LEFT OUTER JOIN v USING (acc) LEFT OUTER JOIN j USING (acc) LEFT OUTER JOIN seq USING (acc);"
    print(query)
    cur.execute(query)

    # *************** Write result of queries to a file ***************

    # Dump info from all_info table to a file
    result = cur.execute('SELECT * FROM all_info')
    header = "\t".join([description[0] for description in result.description])
    colnames = header.split()
    i_v = 0
    i_j = 0
    # Check which column contains the V_gene
    for i in range(len(colnames)):
        if colnames[i] == "V_gene":
            i_v = i
        elif colnames[i] == "J_gene":
            i_j = i
    # Print table to file
    for row in result:
        str_row = list()
        for i in range(len(row)):
            str_row.append(str(row[i]))

        # Make changes to the V_gene name
        str_row[i_v], V_sub, V_main = clean_name(str_row[i_v])

        # Make changes to the J_gene name
        str_row[i_j], J_sub, J_main = clean_name(str_row[i_j])

        # Print the entry
        print("\t".join(str_row) + "\t" + V_sub + "\t" + J_sub + "\t" + V_main, file=fhOut)
    fhOut.close()

    # Import table all_info again and then create clone list
    colnames.append("V_sub")
    colnames.append("J_sub")
    colnames.append("V_main")
    query = "DELETE FROM all_info"
    print(query)
    cur.execute(query)
    query = "ALTER TABLE all_info ADD COLUMN V_sub"
    print(query)
    cur.execute(query)
    query = "ALTER TABLE all_info ADD COLUMN J_sub"
    print(query)
    cur.execute(query)
    query = "ALTER TABLE all_info ADD COLUMN V_main"
    print(query)
    cur.execute(query)
    import_data(args.outfile, "\t", "all_info", colnames)

    # Count number of different V and J genes assigned to one accession code
    query = "CREATE TABLE accs_v_j AS SELECT acc, COUNT(DISTINCT V_main) AS nr_v_mains, COUNT(DISTINCT V_sub) AS nr_v_subs, COUNT(DISTINCT V_gene) AS nr_v_alleles, COUNT(DISTINCT J_sub) AS nr_j_subs, COUNT(DISTINCT J_gene) AS nr_j_alleles FROM all_info GROUP BY acc"
    print(query)
    cur.execute(query)

    # Count number of different sites found in one accession code
    query = "CREATE TABLE sites_count AS SELECT acc,readingframe,COUNT(*) AS nr_sites FROM sites GROUP BY acc,readingframe"
    print(query)
    cur.execute(query)

    # Combine all_info with v and j counts per accession
    query = "CREATE TABLE all_info_nrs AS SELECT all_info.*,accs_v_j.*,sites_count.* FROM all_info JOIN accs_v_j USING (acc) LEFT JOIN sites_count ON all_info.acc=sites_count.acc AND all_info.readingframe=sites_count.readingframe"
    print(query)
    cur.execute(query)

    # Write final all_info table to a file
    fhOut = open(args.outfile, 'w')   # This will overwrite the original all_info file
    result = cur.execute('SELECT * FROM all_info_nrs ORDER BY acc')
    print("\t".join([description[0] for description in result.description]), file=fhOut)
    for row in result:
        str_row = list()
        for i in range(len(row)):
            str_row.append(str(row[i]))
        print("\t".join(str_row), file=fhOut)
    fhOut.close()

    # ************ Make clone reports and write them to a file ****************

    # Create a clone report based on V, J and CDR3peptide
    query = "CREATE TABLE clones AS SELECT V_gene, J_gene, cdr3pep, count(DISTINCT acc) AS freq, count(DISTINCT beforeMID) AS uniq_umis FROM all_info_nrs WHERE V_gene!='None' AND J_gene!='None' and cast(cdr3_qual_min as int)>=" + str(args.qual) + " and (cast(V_flag as int)=0 or cast(V_flag as int)=16) and (cast(J_flag as int)=0 or cast(J_flag as int)=16) GROUP BY V_gene, J_gene, cdr3pep"
    print(query)
    cur.execute(query)

    result = cur.execute("SELECT SUM(freq) AS total_reads, SUM(uniq_umis) AS total_umis FROM clones")
    for row in result:
        if row[0] is None or row[1] is None:
            print("No entries in clones")
            exit()
        total_reads = int(row[0])
        total_umis = int(row[1])

    # Write clones to a file
    result = cur.execute('SELECT * FROM clones ORDER BY freq DESC')
    print("\t".join([description[0] for description in result.description]) + "\tread_perc\tumi_perc", file=fhClones)
    for row in result:
        str_row = list()
        for i in range(len(row)):
            str_row.append(str(row[i]))
        # Calculate percentage
        read_perc = 100 * float(str_row[3]) / float(total_reads)  # column 3 is the frequency
        str_row.append(str(read_perc))
        umi_perc = 100 * float(str_row[4]) / float(total_umis)  # column 3 is the frequency
        str_row.append(str(umi_perc))
        print("\t".join(str_row), file=fhClones)
    fhClones.close()

    # Create a clone report based on V_sub, J_sub and CDR3peptide
    query = "CREATE TABLE clones_subs AS SELECT V_sub, J_sub, cdr3pep, count(DISTINCT acc) AS freq, count(DISTINCT beforeMID) AS uniq_umis FROM all_info_nrs WHERE V_sub!='None' AND J_sub!='None' and cast(cdr3_qual_min as int)>=" + str(args.qual) + " and (cast(V_flag as int)=0 or cast(V_flag as int)=16) and (cast(J_flag as int)=0 or cast(J_flag as int)=16) GROUP BY V_sub, J_sub, cdr3pep"
    print(query)
    cur.execute(query)

    result = cur.execute("SELECT SUM(freq) AS total_reads, SUM(uniq_umis) AS total_umis FROM clones_subs")
    for row in result:
        total_reads = int(row[0])
        total_umis = int(row[1])

    # Write clones to a file
    result = cur.execute('SELECT * FROM clones_subs ORDER BY freq DESC')
    print("\t".join([description[0] for description in result.description]) + "\tread_perc\tumi_perc", file=fhClonesSubs)
    for row in result:
        str_row = list()
        for i in range(len(row)):
            str_row.append(str(row[i]))
        # Calculate percentage
        read_perc = 100 * float(str_row[3]) / float(total_reads)  # column 3 is the frequency
        str_row.append(str(read_perc))
        umi_perc = 100 * float(str_row[4]) / float(total_umis)  # column 3 is the frequency
        str_row.append(str(umi_perc))
        print("\t".join(str_row), file=fhClonesSubs)
    fhClonesSubs.close()

    # Create a clone report based on V_main, J_sub and CDR3peptide
    query = "CREATE TABLE clones_mains AS SELECT V_main, J_sub, cdr3pep, count(DISTINCT acc) AS freq, count(DISTINCT beforeMID) AS uniq_umis FROM all_info_nrs WHERE V_main!='None' AND J_sub!='None' and cast(cdr3_qual_min as int)>=" + str(args.qual) + " and (cast(V_flag as int)=0 or cast(V_flag as int)=16) and (cast(J_flag as int)=0 or cast(J_flag as int)=16) GROUP BY V_main, J_sub, cdr3pep"
    print(query)
    cur.execute(query)

    result = cur.execute("SELECT SUM(freq) AS total_reads, SUM(uniq_umis) AS total_umis FROM clones_mains")
    for row in result:
        total_reads = int(row[0])
        total_umis = int(row[1])

    # Write clones to a file
    result = cur.execute('SELECT * FROM clones_mains ORDER BY freq DESC')
    print("\t".join([description[0] for description in result.description]) + "\tread_perc\tumi_perc", file=fhClonesMains)
    for row in result:
        str_row = list()
        for i in range(len(row)):
            str_row.append(str(row[i]))
        # Calculate percentage
        read_perc = 100 * float(str_row[3]) / float(total_reads)  # column 3 is the frequency
        str_row.append(str(read_perc))
        umi_perc = 100 * float(str_row[4]) / float(total_umis)  # column 3 is the frequency
        str_row.append(str(umi_perc))
        print("\t".join(str_row), file=fhClonesMains)
    fhClonesMains.close()

    # ************ Totals ************

    result = cur.execute('select count(*) from all_info_nrs')
    for row in result:
        print("Total rows in all_info:", row[0], file=fhTotal)

    result = cur.execute('select count(distinct acc) from all_info_nrs')
    for row in result:
        print("Unique reads in all_info_nrs:", row[0], file=fhTotal)

    result = cur.execute("select count(*) from all_info_nrs where V_gene!='None' and J_gene!='None'")
    for row in result:
        print("Total rows with V and J in all_info:", row[0], file=fhTotal)

    result = cur.execute("select count(distinct acc) from all_info_nrs where V_gene!='None' and J_gene!='None'")
    for row in result:
        print("Unique reads with V and J in all_info:", row[0], file=fhTotal)

    result = cur.execute('SELECT SUM(freq) FROM clones')
    for row in result:
        print("Total reads in clones table:", row[0], file=fhTotal)

    result = cur.execute('SELECT SUM(freq) FROM clones_subs')
    for row in result:
        print("Total reads in clones_subs table:", row[0], file=fhTotal)

    result = cur.execute('SELECT SUM(freq) FROM clones_mains')
    for row in result:
        print("Total reads in clones_mains table:", row[0], file=fhTotal)

    fhTotal.close()

    # Close the connection
    con.close()

    print("Finished")

    exit()
