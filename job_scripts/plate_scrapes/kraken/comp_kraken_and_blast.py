

import argparse
import sys

from mypyli import taxstring, kraken

def parse_kraken(kraken_f):
    id_dict = {}
    krak = kraken.KrakenRecord.parse_kraken_file(kraken_f, iterate=True)

    for entry in krak:
        if entry.classified:
            id_dict[entry.name] = entry.tax_id
        else:
            id_dict[entry.name] = "unclassified"

    return id_dict

def parse_blast6(blast6_f):
    id_dict = {}
    with open(blast6_f, 'r') as IN:
        for line in IN:
            elements = line[:-1].split("\t")
            contig_gene = elements[0]
            contig = contig_gene.split("|")[0]

            taxid = elements[-1]
            # take the first assignment if there are more than one
            if ";" in taxid:
                taxid = taxid.split(";")[0]

            if contig in id_dict:
                id_dict[contig][taxid] = id_dict[contig].get(taxid, 0) + 1
            else:
                id_dict[contig] = {taxid: 1}


    # make a final dictionary using the most abundanct taxid for each contig
    final_dict = {}
    for contig, tid_dicts in id_dict.items():

        tids = sorted(tid_dicts, key=lambda tid: tid_dicts[tid], reverse=True)
        final_dict[contig] = tids[0]


    return final_dict

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-k", help="kraken output", required=True)
    parser.add_argument("-b", help="blast output", required=True)

    args = parser.parse_args()

    kraken_dict = parse_kraken(args.k)
    blast_dict = parse_blast6(args.b)


    both_unclass = 0
    k_unclass = 0
    b_unclass = 0
    same_class = 0
    lin_divides = {'k': 0, 'p': 0, 'c': 0, 'o': 0, 'f':0, 'g': 0, 's': 0, 'N': 0}
    
    correct_contigs = []
    processed = 0
    total = len(kraken_dict)
    for contig, k_assign in kraken_dict.items():
        b_assign = blast_dict.get(contig, "unclassified")

        if k_assign == "unclassified":
            if b_assign == "unclassified":
                both_unclass += 1
            else:
                k_unclass += 1
        else:
            if b_assign == "unclassified":
                b_unclass += 1
            else:
                try:
                    k_tax = taxstring.TaxString(tax=k_assign, is_id=True, lookup=True)
                    b_tax = taxstring.TaxString(tax=b_assign, is_id=True, lookup=True)

                    if k_tax.get_tax_string() == b_tax.get_tax_string():
                        same_class += 1
                        correct_contigs.append(contig)
                    else:
                        divide = taxstring.TaxString.lineage_divides_at(k_tax, b_tax)
                        lin_divides[divide[0]] += 1
                        if divide in ["family", "genus", "species"]:
                            correct_contigs.append(contig)
                except: # no acceptions allowed
                    print("Exception for contig: {}".format(contig))

        if processed / 500 == int(processed / 500):
            print("{} of {}({}%) completed...".format(processed, total, int(processed / total * 100)))
        processed += 1
        #if stop > 10:
        #    break
        #stop += 1

    print("\n"*3)
    print(("both_unclass", both_unclass))
    print(("k_unclass", k_unclass))
    print(("b_unclass", b_unclass))
    print(("same_class", same_class))

    print("Splits at:")
    for k in ['k', 'p', 'c', 'o', 'f', 'g', 's', 'N']:
        print("  {}\t{}".format(k, lin_divides[k]))

    with open("correct_contigs.txt", 'w') as OUT:
        for contig in sorted(correct_contigs):
            OUT.write("{}\n".format(contig))
