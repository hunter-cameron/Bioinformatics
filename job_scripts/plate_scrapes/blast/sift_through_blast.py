

from mypyli import blastparser
import sys
import os

MIN_ID = 90
MIN_LEN = 1000

for blast_f in sys.argv[1:]:
  
    counts = {}
    with open(blast_f, 'r') as IN:
        for record in blastparser.parse(IN, "6"):
            if record.get_perc_id() >= MIN_ID and record.get_length() >= MIN_LEN:
                subject = record.get_subj()
                isolate = subject.split("_")[0]
            
                counts[isolate] = counts.get(isolate, 0) + 1
            else:
                counts["matchless"] = counts.get("matchless", 0) + 1
    
    print("\t".join([blast_f.replace(".blast6out", "")] + ["{}:{}".format(isolate, str(counts[isolate])) for isolate in counts]))
