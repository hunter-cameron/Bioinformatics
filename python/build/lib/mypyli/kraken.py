
import sys

class KrakenRecord(object):
    
    @classmethod
    def parse_kraken_file(cls, kraken_f, iterate=False):
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
        self.id_hits = [line_dict['id_hits'].split(" ")]
        
    def __repr__(self):
        return "{} classified as TID:{}".format(self.name, self.taxid) if self.classified else "{} not classified".format(self.name)




if __name__ == "__main__":

    records = KrakenRecord.parse_kraken_file(sys.argv[1])
    [print(record) for record in records]

    KrakenRecord.print_count_matrix(records)
