

from Bio import SeqIO
import sys
import re

def get_entries_by_locid(gbk_file, id_list):
    print(id_list)
    results = dict()
    for record in SeqIO.parse(open(gbk_file,"r"), "genbank"):       # record == contig
        for feature in record.features:                             # all annotations per contig
            if feature.type == "CDS":                               # get only CDS
                #sys.exit()
                if "locus_tag" in feature.qualifiers:               # check if CDS has a locus tag (it should)
                    if feature.qualifiers['locus_tag'][0] in id_list:  # check if locus tag is on the list
                        results[feature.qualifiers['locus_tag'][0]] = {
                            "location": feature.location,
                            "product": feature.qualifiers['product'][0],

                        }
                #sys.exit()
    return results
        

def read_locid_list(id_file):
    """ Returns a list of sorted ids from a file """
    with open(id_file, 'r') as IN:
        return sorted([line[:-1] for line in IN])





if __name__ == "__main__":
   
    id_file = sys.argv[1]
    gbk_file = sys.argv[2]
  
    id_list = []

    with open(id_file, 'r') as IN:
        for line in IN:
            qry = line.split("\t")[0]
            loctag = qry.split(" ")[0]
            id_list.append(loctag)
    
    
    id_info = get_entries_by_locid(gbk_file, id_list)

    for tag in id_list:
        if tag not in id_info:
            print("Locus tag '{}' not found.".format(tag))


    
    with open(id_file, 'r') as IN, open("final_matrix.txt", 'w') as OUT:
        for line in IN:
            if line.startswith("qry"):
                OUT.write("\t".join(["locus_tag", "contig", "contig_length", "start", "end", "strand", "product", "closest_match", "perc_id", "aln_length", "query_cov", "closest_match_cov", "bitscore"]) + "\n")

            else:
                elements = line[:-1].split("\t")

                qry_info = elements[0].split(" ")
                locid = qry_info[0]
                contig = qry_info[2]
                m = re.search("len_(?P<length>\d+)_", contig)
                contig_length = m.group("length")


                OUT.write("\t".join([locid, contig, contig_length, str(id_info[locid]['location'].start), str(id_info[locid]['location'].end), str(id_info[locid]['location'].strand), id_info[locid]['product'], "\t".join(elements[1:])]) + "\n")

