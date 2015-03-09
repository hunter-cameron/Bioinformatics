

from mypyli import samparser, kraken, taxstring

import argparse
import pickle
import sys

def get_tax_assignments(kraken_f):
    """
    Parse the kraken output file and return each read's assignment.
    Often only one of the pairs gets assigned. Make sure to take the one that gets assigned.
    If both get assigned to different taxa, we will just keep the first.
    """
    
    tax_dict = {}
    for krakenr in kraken.KrakenRecord.parse_kraken_file(kraken_f, iterate=True):
        #print(str(krakenr.classified))
        if tax_dict.get(krakenr.name, 0):
            if tax_dict[krakenr.name] != "unclassified":
                continue

        if krakenr.classified:
            tax_dict[krakenr.name] = krakenr.tax_id
        else:
            tax_dict[krakenr.name] = "unclassified"
    return tax_dict

def get_read_alignments(sam_f):
    """
    Reads a SAM file and returns a data structure that links the contigs to reads.
    Alternatively, could write it to temp file to keep it out of memory until it is needed and return the path.
    """
    sparser = samparser.SamParser(sam_f=sam_f, aligned_only=True, mapq=20, mismatches=1)
    
    # parse all the hits into this to make sure multi mapping hits map to the same contig
    hit_dict = {}
    ambig_reads = 0
    processed_reads = 0
    for hit in sparser.parse_sam_file():
        processed_reads += 1
        if hit_dict.get(hit['qname'], 0):
            if hit_dict[hit['qname']] != hit['rname']:
                print("Warning read: {} aligns to two different contigs".format(hit['qname']), file=sys.stderr)
                ambig_reads += 1
            else:
                continue
        else:
            hit_dict[hit['qname']] = hit['rname']

    print("{} of {} processed reads were ambiguous.".format(ambig_reads, processed_reads))

    # condense the hit dict into a contig dict
    contig_dict = {}
    for read, contig in hit_dict.items():
        if contig_dict.get(contig, 0):
            contig_dict[contig].append(read)
        else:
            contig_dict[contig] = [read]

    return contig_dict

def link_contigs_with_tax(contig_dict, tax_dict):
    linked_dict = {}
    for contig, read_arr in contig_dict.items():
        linked_dict[contig] = {}
        for read in read_arr:
            tax = tax_dict.get(read)
            if tax:
                linked_dict[contig][tax] = linked_dict[contig].get(tax, 0) + 1

            else:
                print("Read {} not found in kraken output".format(read), file=sys.stderr)

    return linked_dict

def lookup_tax_by_id(linked_dict):
    """
    Looks up a tax ids found in a dict of structure: linked_dict = {contig: {taxid: count}}
    Returns a dict = taxid: taxonomy
    """

    #tax_ids = {k: 0 for k in [taxa for taxa in linked_dict.values()]}

    tax_ids = {}
    for taxa in linked_dict.values():
        for tid in taxa:
            tax_ids[tid] = 0

    #tax_values = {}
    for taxid in tax_ids:
        if taxid == "unclassified":
            tax_ids[taxid] = "unclassified"
        else:
            print(taxid)
            taxstr = taxstring.TaxString(tax=taxid, is_id=True, lookup=True)
            tax_ids[taxid] = taxstr.get_tax_string()
     #       if tax_values.get(taxstr.get_tax_string(), 0):
      #          print(("Same", tax_values[taxstr.get_tax_string()], taxid))
                #sys.exit()
       #     else:
        #        tax_values[taxstr.get_tax_string()] = [taxid]

    return tax_ids

def print_kraken_otu_table(linked_dict, tax_ids, out="contig_taxonomy_table.txt"):
    """
    Prints an OTU-like table of reads and taxonomies assigned to contigs.
    Because taxon ids are not necessarily unique, must check and merge same taxonomies.
    """
    unique_tax = get_unique_dict(tax_ids)


    sorted_tax =  sorted(unique_tax.keys())
    with open(out, 'w') as OUT:
        OUT.write("\t".join(["contig"] + sorted_tax) + "\n")

        
        for contig, tid_dict in linked_dict.items():
            to_print = []
            to_print.append(contig)

            
            for tax in sorted_tax:
                count = 0
                for tid in unique_tax[tax]:
                    count += tid_dict.get(tid, 0)

                to_print.append(str(count))

            OUT.write("\t".join(to_print) + "\n")

def get_unique_dict(tax_dict):
    
    print("Eliminating duplicate entries...", file=sys.stderr)
    unique_dict = {}
    for tid, tax in tax_dict.items():
        unique_dict[tax] = unique_dict.get(tax, []) + [tid]
        
    return unique_dict

        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Links a Kraken assignment file with a sam file to assign phylogeny to contigs")
    parser.add_argument("-kraken", help="output file from kraken")
    parser.add_argument("-sam", help="sam file from aligning the reads to the contigs")
    parser.add_argument("-out", help="the output path for the pickled data structure", default="kraken2contig_pickle.txt")
    parser.add_argument("-mode", help="the mode in which to run. The translate option requires the ability to conect to the NCBI servers (i.e. may not work with bsub).", required=True, choices=["link", "print", "both"], default="link")
    parser.add_argument("-k2c", help="the kraken to contig pickle file")
    args = parser.parse_args()

    # determine the mode and check for required arguments
    if args.mode in ["link", "both"]:
        if not args.kraken or not args.sam:
            parser.error("selected -mode option requires options -kraken and -sam")

        contig_dict = get_read_alignments(args.sam)
        #print(contig_dict)

        tax_dict = get_tax_assignments(args.kraken)
        #print(tax_dict)

        linked_dict = link_contigs_with_tax(contig_dict, tax_dict)
   
        # dump the pickle either way to avoid having to run it all again
        with open(args.out, 'wb') as OUT:
            pickle.dump(linked_dict, OUT)

        if args.mode == "both":
            # try to free up some memory before going on to the next step
            contig_dict = None
            tax_dict = None


    if args.mode in ["print", "both"]:
        if args.mode == "print":
            if not args.k2c:
                parser.error("selected -mode option requires option -k2c")

            else:
                linked_dict = pickle.load(open(args.k2c, 'rb'))


        tax_ids = lookup_tax_by_id(linked_dict)
        print_kraken_otu_table(linked_dict, tax_ids)
