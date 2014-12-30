
import sys

class KrakenIO(object):
    """ This is basically just a class to mimic SeqIO from Biopython to give user friendly access """ 
    
    @classmethod
    def parse(cls, kraken_fh):
        """
        Class method that generates KrakenRecord Objects. 
        
        """

        record_list = []
        for line in kraken_fh:
            yield KrakenRecord(line)


    def __init__(self):
        pass



class KrakenRecord(object):
    
    tree = TaxTree.load_tree("/nas02/home/h/j/hjcamero/scripts/python/mypyli/taxtree_kingdom.pickle")

    @classmethod
    def parse_kraken_file(cls, kraken_f, iterate=True):
        """
        Class method that returns a list of all kraken entries in the output 
        
        Also has the option to work as a generator and iterate back out objects (to selectively keep some)
        """

        record_list = []
        with open(kraken_f) as IN:
            if iterate:
                for line in IN:
                    yield cls(line)
            else:
                return [cls(line) for line in IN]

    @staticmethod
    def print_count_matrix(record_list, outfile="kraken_taxonomy_counts.txt"):
        counts = {}
        for record in record_list:
            if record.classified:
                counts[record.tax_id] = counts.get(record.tax_id, {'count':0, 'length':0})
                counts[record.tax_id]['count'] += 1
                counts[record.tax_id]['length'] += int(record.length)
            else:
                counts['unclassified'] = counts.get('unclassified', {'count':0, 'length':0})
                counts['unclassified']['count'] += 1
                counts['unclassified']['length'] += int(record.length)

        with open(outfile, 'w') as OUT:
            print("taxid\tcount\tavg_length")
            for key, value  in counts.items():
                OUT.write("\t".join([str(key), str(value['count']), str(value['length'] / value['count'])]) + "\n")


    def __init__(self, line):
        # "chomp" the line
        line = line.rstrip()
        line_dict = {key: value for key, value in zip(["classified", "name", "tax_id", "length", "id_hits"], line.split("\t"))}

        # set classified to true if C, false otherwise
        self.classified = line_dict['classified'] == 'C'
        self.name = line_dict['name']
        self.taxid = line_dict['tax_id']
        self.length = line_dict['length']
        self.kmer_hits = line_dict['id_hits']   # leave in raw form in order for faster processing when user doesn't have about hits
        
    def __repr__(self):
        return "{} classified as TID:{}".format(self.name, self.taxid) if self.classified else "{} not classified".format(self.name)


    def get_kraken_confidence(self):
        """ This is a confidence measure that kraken talks about in their readme.

        It is simply the # of kmers mapping to nodes rooted in the assigned node divided by
        the number of non-ambiguous sequences (assigned to A:)
        """
        if not self.classsified:
            return 0
        
        kmer_dict = self._convert_hits_to_dict()

        accurate_count = 0
        inaccurate_count = 0
        tree = KrakenRecord.tree
        assigned = tree.lookup_taxid(self.taxid)
        for taxid in kmer_dict:
            current = tree.lookup_taxid(taxid)

            if current == 'A' or current == "0":
                continue

            elif assigned == current or assigned.is_parent_of(current):
                accurate_count += int(kmer_dict[taxid])
            else:
                inaccurate_count += int(kmer_dict[taxid])
                

        

    def _convert_hits_to_dict(self):
        """ This converts the hits to a dict. But I'm leaving options open to convery to some sort of list because I might want to preserve the order in which kmers mapped to find misassembled contigs. """
        if type(self.kmer_hits) is str:
            kmer_dict = {}
            for hit in kmer_hits.split(" "):
                [taxid, count] = hit.split(":")
                kmer_dict[taxid] = kmer_dict.get(taxid, 0) + int(count)

            return kmer_dict

if __name__ == "__main__":

    records = KrakenRecord.parse_kraken_file(sys.argv[1])
    [print(record) for record in records]

    KrakenRecord.print_count_matrix(records)
