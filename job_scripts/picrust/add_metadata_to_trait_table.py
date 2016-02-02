
import argparse
import picrust

def load_metadata(metadata_f):
    metadata = {}
    with open(metadata_f, 'r') as IN:
        headers = IN.readline().rstrip().split("\t")[1:]

        for line in IN:
            fields = line.rstrip().split("\t")

            for index, field in enumerate(fields[1:]):
                try:
                    # store in metadata at trait -> metadata header
                    metadata[headers[index]][fields[0]] = field
                except KeyError:
                    metadata[headers[index]] = {fields[0]: field}

    return metadata

def add_metadata(trait_table_f, metadata, out):
    ttm = picrust.TraitTableManager(trait_table_f)

    traits = ttm.get_ordered_traits()
    with open(out, 'w') as OUT:
        OUT.write("\t".join([ttm.entry_header] + traits) + "\n")

        for entry in ttm:
            ttm.write_entry(entry, OUT, traits)

        for name in metadata:
            tte = picrust.TraitTableEntry(name)
            tte.traits = metadata[name]

            ttm.write_entry(tte, OUT, traits)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Adds the metadata from the precalculated files to a trait table to enable the usage of catagorize_by_function.py PICRUSt script.")
    parser.add_argument("-m", help="trait metadata file from deconstructed PICRUSt", required=True)
    parser.add_argument("-t", help="trait table to add the metadata to", required=True)
    parser.add_argument("-o", help="path for new trait table", required=True)

    args = parser.parse_args()

    metadata = load_metadata(args.m)
    add_metadata(args.t, metadata, args.o)

