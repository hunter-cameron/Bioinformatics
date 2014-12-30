#!/usr/bin/env python3

from Bio import SeqIO
import argparse
import os
import sys

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("path", nargs="+", help="Paths to the files/directories (default=file).")
    parser.add_argument("-d","--dir", help="Treat paths as directories.", action="store_true")

    args = parser.parse_args()

    #check if directories are directories
    if args.dir == True:
        for path in args.path:
            if (not os.path.isdir(path)):
                sys.exit("Path: {} is not a valid directory".format(path))
            else:
                files = get_files(path)
                #print(files)
                gbk2faa(files)

    #check if files are files
    else:
        for path in args.path:
            if (not os.path.exists(path)):
                sys.exit("Path: {} is not a valid path.".format(path))

        gbk2faa(args.path)

def gbk2faa(files):
    for file in files:

        try:

            out_file = get_filename(file) + ".faa"
            print("Writing to: {}".format(out_file))
            OUT = open(out_file, "w")
            for record in SeqIO.parse(file, "genbank"):
                for feature in record.features:
                    if feature.type == "CDS":
                        assert len(feature.qualifiers['translation']) == 1
                        OUT.write(">{} from {} product={}\n{}\n".format(feature.qualifiers['locus_tag'][0], record.name, feature.qualifiers['product'][0], feature.qualifiers['translation'][0]))
            OUT.close()
        
        except Exception as e:
            print("Could not process {}".format(file))
            print(e)




def get_files(dir):
    gbk = []
    for name in os.listdir(dir):
        #print(os.path.splitext(name))
        if os.path.splitext(name)[1] in [".gbk", ".GBK"]:
            gbk.append(dir + name)
    return gbk

def get_filename(path):
    """
    Function to roll together basename and splitext to remove the path prefix and the extension
    """

    name_ext = os.path.basename(path)
    name = os.path.splitext(name_ext)[0]
    return name





if __name__ == "__main__":
    main()
