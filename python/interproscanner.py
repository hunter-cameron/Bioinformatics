#!/usr/bin/env python3

class OutputParser(object):

    def __init__(self, fmt, file_to_parse):
        self.fmt = fmt.upper()
        self.file_to_parse = file_to_parse

    def parse(self):
        """ Wrapper function in case I want to be able to parse other formats at some time """

        if self.fmt == "GFF3":
            return self.parse_gff3()


    def parse_gff3(self):
        """ Function to handle parsing of the GFF3 file format, returns a list of Query objects """

        queries = {}
        with open(self.file_to_parse, 'r') as IN:
            for line in IN:
                if line.startswith("##"):
                    if line.startswith("##FASTA"):
                        break
                    continue

                #elements = line[:-1].split("\t")
                # make a dictionary of the elements in the line
                hit_dict = {key: value for key, value in zip(["query", "source", "type", "start", "end", "score", "strand", "phase", "attributes"], line[:-1].split("\t"))}

                if hit_dict["type"] == "polypeptide":
                    queries[hit_dict["query"]] = Query(hit_dict)
                if hit_dict["score"] == ".":
                    continue

                if hit_dict["query"] in queries:
                    queries[hit_dict["query"]].add_hit(hit_dict)
                else:
                    print("WARNING! Not found in dict")
                    sys.exit()
        
        return queries.values()

class Query(object):

    def __init__(self, hit_dict):

        self.length = int(hit_dict["end"]) - int(hit_dict["start"])
        self.name = hit_dict["query"]
        self.hits = []

    def add_hit(self, hit_dict):
        """ Need to move most of this processing to the GFF stuff """

        elements = hit_dict["attributes"].split(";")
        
        ## Need to make sure all the params I want here are defined so I don't run into issues later
        hit_dict["go_terms"] = []

        for element in elements:
            if element.startswith("Name="):
                hit_dict["subject"] = element[5:]
            if element.startswith("signature_desc="):
                hit_dict["desc"] = element[15:]
            if element.startswith("Ontology_term="):
                element = element.replace("\"", "")
                terms = element[14].split(",")
                hit_dict["go_terms"] = element[14:].split(",")

    
        # convert all number like things to numbers
        for key in hit_dict:
            if not isinstance(hit_dict[key], list):
                try:
                    hit_dict[key] = float(hit_dict[key])
                except ValueError:
                    continue

        self.hits.append(hit_dict)

    def remove_bad_hits(self, max_e=.05, min_cov=.5):
        hits_to_keep = []
        for hit in self.hits:
            if hit["score"] > max_e:
                continue
            
            if (hit["end"] - hit["start"]) < min_cov * self.length:
                continue

            hits_to_keep.append(hit)

        self.hits = hits_to_keep

    def get_all_GO_terms(self):
        """ Returns a list with all the GO terms for the hits. Do this after remove_bad_hits """

        terms = {}
        for hit in self.hits:
            terms[hit["go_terms"]] == 1

        return list(terms.keys())

    def get_all_pfams(self):
        """ Returns a list with all the pfams the query has hits to """

        terms = {}
        for hit in self.hits:
            if hit["source"] == "Pfam":
                terms[hit["subject"]] = 1

        return list(terms.keys())

class InterproScanner(object):

    def __init__(self,
            # interproscan vars
            fasta_in,
            out_path="interproscan_results.txt",
            fmt="GFF3",
            bin="/proj/dangl_lab/apps/my_interproscan/interproscan-5.3-46.0/interproscan.sh",
            othr_args="",

            # bsub vars 
            threads=8,
            queue="week",
            stdout="interproscanner_bsub.out",
            stderr="interproscanner_bsub.err",

            # batch vars
            seqs_per_file=25
            ):


        self.ips_vars = {
                'fasta_in': [fasta_in],
                'out_path': out_path,
                'fmt': fmt,
                'bin': bin,
                'other_args': othr_args
                }

        self.bsub_vars = {
                'threads': threads,
                'queue': queue,
                'stdout': stdout,
                'stderr': stderr
                }

        self.batch_vars = {
                'seqs_per_file'
                }
        

    def run_interproscan(self):
        

        # count sequences
        seqs = 0
        with open(self.ips_vars["fasta_in"][0], 'r') as IN:
            for line in IN:
                if line.startswith(">"):
                    seqs += 1
        
        # split fasta if needed
        fas_files = []
        if seqs > self.batch_vars["seqs_per_file"]:
            fas_files = self._split_fasta()

        # run command on each fasta
        for fas_file in fas_files:
            command = self._build_command(fas_file)
            print("Executing: " + command)


    def _split_fasta(self):
        """ Splits a the fasta into multiple parts and changes the fasta in to the parts """

        counter = 1
        file_indx = 1
        fastas = []

        with open(self.ips_vars["fasta_in"][0], 'r') as IN:
            OUT = open("tmp_iproscanner_{}.fasta".format(file_indx), 'w')
            fastas.append("tmp_iproscanner_{}.fasta".format(file_indx))

            for line in IN:
                # increment counter for each header
                if line.startswith(">"):
                    counter += 1
                    
                    # reset and open new file if counter is enough
                    if counter > self.batch_vars["seqs_per_file"]:
                        counter = 1
                        file_indx += 1

                        OUT.close()
                        OUT = open("tmp_iproscanner_{}.fasta".format(file_indx), 'w')
                        fastas.append("tmp_iporscanner_{}.fasta".format(file_indx))
                        
                OUT.write(line)

        self.fasta_in = fastas
            
    def _build_command(self, fasta_file):
        """ Builds a command to run interproscan for a given fasta file """

        # shell with out and err
        command = "bsub -o {} -e {}".format(self.bsub_params["stdout"], self.bsub_params["stderr"])

        # add threading
        hosts = self.bsub_params["threads"] // 8    # use as few hosts as possible
        command += "-n {} -R 'span[hosts={}']".format(self.bsub_params["threads"], hosts)

        # add interpro with base options and file
        command += "{} --goterms -dp -i {}".format(self.ips_params["bin"], fasta_file)

        # add output options
        command += "-o {} -f {}".format(self.ips_params["out_path"]. self.ips_params["fmt"])

        # add any other options
        command += self.ips_params["other_args"]

        return command

def go_term_enrichment(queries):

    term_counts = {}
    for qry in queries:

        qry.remove_bad_hits(max_e=.05, min_cov=.01)
        terms = qry.get_all_GO_terms()
        
        for term in terms:
            term_counts[term] = term_counts.get(term, 0) + 1


    for term in sorted(term_counts, key=lambda key: term_counts[key], reverse=True):
        print(term + "\t" + str(term_counts[term]))

def pfam_dump(queries):
    for qry in queries:

        qry.remove_bad_hits(max_e=.05, min_cov=.01)
        #print("Getting pfams for " + qry.name)
        pfams = qry.get_all_pfams()

        for pfam in pfams:
            print(pfam + "\t" +  qry.name)


if __name__ == "__main__":

    import sys

    my_file = sys.argv[1]

    outp = OutputParser(fmt="GFF3", file_to_parse=my_file)

    queries = outp.parse()

    pfam_dump(queries)
