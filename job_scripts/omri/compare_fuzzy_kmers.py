
import argparse
import os
import sys
import math
import time
from mypyli import jellyfish
import numpy as np
import logging
import re
import pickle
from Bio import SeqIO
import queue
import multiprocessing
import psutil
import sqlite3
#import matplotlib.pyplot as plt

import objgraph
from pympler import tracker, muppy, summary

#from matplotlib import rcParams
#rcParams.update({'figure.autolayout': True})


logging.basicConfig()
LOG = logging.getLogger(__name__)
LOG.setLevel("DEBUG")

"""

TODO: Rewrite using sqlite. Should be faster and use less memory. Basically, I have written a poor version of SQLite using my dumping method.

TODO: The mergesort algorithm may be incorrect. For amplicons, it reports a correct number of locations but has more amplicons (not too many more)

    For Kmers, it reports less than half as many results (locations) though there were a consistent number (4911) kmers conserved between all genomes.

TODO: Print some final stats about amplicons (number unique, etc).

TODO: It would be great to be able to be able to load the indexes from a previous run to avoid having to recalculate them.
"""

class KmerCounterSQL(object):


    def __init__(self, database, k, threads, mismatches=1, stable_bp=2):
        self.database_f = database
        self.k = k
        self.threads = threads
        self.mismatches = mismatches

        # stables are # bases at beginning and end that cannot be mismatches
        self.stable_bp = stable_bp

        # open/make the database
        if os.path.isfile(self.database_f): 
            self.open_database()
        else:
            self.open_database()
            self.make_SQL_tables()

    def open_database(self):
        self.database = sqlite3.connect(self.database_f, check_same_thread=False)
        self.dbc = self.database.cursor()

    def make_SQL_tables(self):
        """ Makes SQL tables corresponding to all the info I want """

        # create table that holds the genome names
        self.dbc.execute("""CREATE TABLE Genomes (
                        id INTEGER PRIMARY KEY NOT NULL,
                        name TEXT
                        );
                        """)

        # create table that holds contigs
        self.dbc.execute("""CREATE TABLE Contigs (
                        id INTEGER PRIMARY KEY NOT NULL,
                        name TEXT,
                        genome INTEGER NOT NULL,
                        FOREIGN KEY(genome) REFERENCES Genomes(id)
                        );
                        """)

        # create table that holds all the unique Kmers that are stored in locations
        # id is an integer that results from converting ACGT to base4 and then converting to base10
        self.dbc.execute("""CREATE TABLE Kmers (
                        id INTEGER PRIMARY KEY ON CONFLICT IGNORE NOT NULL,
                        name TEXT
                        );
                        """)

        # create table to hold each location (could add a genome field to make lookups ever easier)
        self.dbc.execute("""CREATE TABLE Locations (
                        id INTEGER PRIMARY KEY NOT NULL,
                        contig INTEGER NOT NULL,
                        start INTEGER,
                        strand INTEGER,
                        kmer INTEGER NOT NULL,
                        FOREIGN KEY(contig) REFERENCES Contigs(id)
                        FOREIGN KEY(kmer) REFERENCES Kmers(id)
                        );
                        """)

        #
        ## Create some tables to aid in calculation of conserved kmers
        ## These tables will be populated after kmer counting
        #

        # create a table of fuzzy kmers 
        self.dbc.execute("""CREATE TABLE FuzzyKmers (
                        id INTEGER PRIMARY KEY ON CONFLICT IGNORE NOT NULL,
                        name TEXT
                        );
                        """)

        # create a many to many table that relates each kmer with all the fuzzy kmers it could be a part of
        self.dbc.execute("""CREATE TABLE KmerToFuzzy (
                        id INTEGER PRIMARY KEY NOT NULL,
                        kmer INT,
                        fuzzy INT,
                        UNIQUE(kmer, fuzzy) ON CONFLICT IGNORE
                        FOREIGN KEY(kmer) REFERENCES Kmers(id)
                        FOREIGN KEY(fuzzy) REFERENCES FuzzyKmers(id)
                        );
                        """)

    @classmethod
    def _recurse_all_kmers(cls, k, curseq=""):
        """ Generates all kmers to avoid having to store a large array """
        if len(curseq) == k:
            yield curseq
        else:
            for nt in ["A", "C", "G", "T"]:
                for kmer in cls._recurse_all_kmers(k, curseq + nt):
                    yield kmer
 

    def count_kmers(self, fastas):
        """
        The counter here should do nothing but update the database.

        In multiple other threads the counter should be running counting kmers and should return multiple lists of tuples. The first contains all the unique kmers that need to be in kmers and the second contains Locations for the counted segment.
        """
        input_queue = multiprocessing.Queue()
        results_queue = multiprocessing.Queue(15)

        # spawn workers
        workers = []
        for indx in range(self.threads):
            worker = KCounter("Counter{}".format(indx), self.k, input_queue, results_queue)
            worker.start()
            workers.append(worker)

        for fasta in fastas:
            genome_name = os.path.splitext(os.path.basename(fasta))[0]

            # check if genome is in database abd skip if so
            if self.genome_in_database(genome_name):
                LOG.info("Genome {} already found in database. Skipping.".format(genome_name))
                continue
            else:
                LOG.info("Counting kmers for {}".format(genome_name))

            # add fasta to the database and get the genome index
            genome_id = self._add_genome(genome_name)

            with open(fasta, 'r') as IN:
                # keep track of how much goes in the pipeline so I'll know how much to expect out
                num_batches = 0
                for record in SeqIO.parse(IN, "fasta"):


                    # add contig to the database and get the contig index
                    contig_id = self._add_contig(record.description, genome_id)

                    # split the sequence up into pieces and put the pieces into the input queue
                    for indx in range(0, len(record.seq), 500000):
                        input_queue.put((genome_id, contig_id, indx, str(record.seq[indx:indx+500000])))
                        num_batches += 1


            # get result batches 
            for bat_num in range(num_batches):
                LOG.debug("Getting batch {} of {}".format(bat_num+1, num_batches))
                kmers, locations = results_queue.get()
        
                LOG.debug("Adding kmers to database.")
                self._add_kmers(kmers)
                LOG.debug("Adding locations to database.")
                self._add_locations(locations)
            self.database.commit()
            print(self.database.total_changes)

        for _ in workers:
            input_queue.put("STOP")

        for worker in workers:
            LOG.debug("Waiting for {} to join...".format(worker.name))
            worker.join()

    #
    ## Populating fuzzy kmer table
    #
    def generate_fuzzy_kmers(self):
        """ Populates the FuzzyKmer and KmerToFuzzy tables 
        
        All kmers without entries in the fuzzy tables are queries and results are written to file
        Then, file is iterated over in multiprocessing and batches of fuzzy entries are added to table.

        I use this schema because SQL complains (and sometimes Segfaults) if any sqlite object (including a cursor) is used in a thread other than the originating thread)
        """

        # make a pool to process results
        pool = multiprocessing.Pool(self.threads)

        # select kmers that do not have an entry in the KmersToFuzzy table (ones that have not yet been processed)
        LOG.info("Selecting kmers to fuzzify...")
        cursor = self.dbc.execute("SELECT Kmers.id, Kmers.name from Kmers LEFT JOIN KmerToFuzzy ON kmers.id = KmerToFuzzy.kmer WHERE KmerToFuzzy.kmer is NULL")


        # write results to a temporary file
        # TODO Actually use the tmpfile module to generate this file
        LOG.info("Writting kmers to a temporary file...")
        with open("tmpkmer.dmp", "w") as OUT:
            self._sql_results_to_file(OUT, cursor)

        iter_fuzzy_kmers = pool.imap(self._generate_fuzzy_entries, self._iter_from_file("tmpkmer.dmp", 100000))

        LOG.info("Fuzzifying and inserting fuzzy kmers into database...")
        for fuzzy_kmers, links in iter_fuzzy_kmers:
            LOG.debug("Inserting {} elements into FuzzyKmers...".format(len(fuzzy_kmers)))
            self.dbc.executemany("INSERT OR IGNORE INTO FuzzyKmers (id, name) VALUES (?, ?)", fuzzy_kmers)
            LOG.debug("Inserting {} elements into KmerToFuzzy...".format(len(links)))
            self.dbc.executemany("INSERT INTO KmerToFuzzy (kmer, fuzzy) VALUES (?, ?)", links)

        self.database.commit()
        
    def _generate_fuzzy_entries(self, kmers):
        """ Takes an iterable of kmers (id, name) and returns a list of fuzzy kmers ready to be inserted into the FuzzyKmer table and a list of links between FuzzyKmers and real kmers"""

        fuzzy_kmers = []
        links = []
        for kmer_id, kmer in kmers:
            for fuzzy_kmer in self._fuzzify_kmer(kmer):
                fuzzy_id_str = fuzzy_kmer.replace("A", "0").replace("C", "1").replace("G", "2").replace("T", "3").replace("N", "4")
                fuzzy_id = int(fuzzy_id_str, 5)
            
                fuzzy_kmers.append((fuzzy_id, fuzzy_kmer))
                links.append((kmer_id, fuzzy_id))

        return fuzzy_kmers, links

    def _fuzzify_kmer(self, kmer, current_seq="", mismatches=0):
        """ 
        Generates all fuzzy kmers with a given number of mismatches (by self.mismatches)
        will not generate fuzzy kmers with < self.mismatches because clusters with mismatch < N are included

        However, for best primer finding, it may be useful (though computationally redundant) to have the most error free fuzzy representative. Therefore, I would recommend finding fuzzy kmers adding a progressive amounts of N's as appropriate results aren't found.

        This is a slow point of my program right now. I need to speed up.  
        """
        # check if we are to the end
        if len(current_seq) >= self.k - self.stable_bp:
            # only yield is the appropriate number of mismatches have been allowed
            if mismatches == self.mismatches:
                yield "".join((current_seq, kmer))

            return
        
        # skip if we are still in the stable start
        elif len(current_seq) < self.stable_bp:
            pass

        else:
            
            # check if the current base can be fuzzy (not too many previous errors)
            if mismatches < self.mismatches:
                # run the ambiguous result
                for result in self._fuzzify_kmer(kmer[1:], current_seq+"N", mismatches+1):
                    yield result
 
            # base cannot be fuzzy, already enough errors
            else:
                # once again, not sure why I can't just return this...
                yield "".join((current_seq, kmer))
                return
        
        # run the next iteration if not returned
        for result in self._fuzzify_kmer(kmer[1:], current_seq+kmer[0], mismatches):
            yield result


    #
    ## Database Commands
    #

    def genome_in_database(self, genome_name):
        """ Checks if a genome is present in the Genomes table """
        cursor = self.dbc.execute("SELECT 1 FROM Genomes WHERE name = '{}' LIMIT 1".format(genome_name))
        results = cursor.fetchone()
        print(results)
        if results:
            return True
        else:
            return False

    def _change_sync(self, value):
        """ Switches pragma sync to either OFF or NORMAL """
        if value in ["OFF", "NORMAL"]:
            self.dbc.execute("PRAGMA synchronous = {}".format(value))
        else:
            raise ValueError("{} is not an appropriate sync mode.".format(value))

    def _add_genome(self, genome_name):
        """ Adds a genome name and returns genome id in database """

        self.dbc.execute("INSERT INTO Genomes (name) VALUES (?)", (genome_name,))
        return self.dbc.lastrowid

    def _add_contig(self, contig_name, genome_id):
        """ Adds a contig and returns the contig id in the database """

        self.dbc.execute("INSERT INTO Contigs (name, genome) VALUES (?, ?)", (contig_name, genome_id))
        return self.dbc.lastrowid

    def _add_kmers(self, kmers):
        """ Insert ignore new kmers from a list of kmer iterables """
        self.dbc.executemany('INSERT OR IGNORE INTO Kmers (id, name) VALUES (?, ?)', kmers)

    def _add_locations(self, locations):
        """ Insert new locations from a list of location tuples """
        self.dbc.executemany('INSERT INTO Locations (contig, start, strand, kmer) VALUES (?, ?, ?, ?)', locations)

    #
    ## Querying the database 
    #
    
    @staticmethod
    def _sql_results_to_file(fh, cursor):
        """ Writes results to a tab delimited file """
        while True:
            results = cursor.fetchmany(10000)
            if not results:
                break
            else:
                for result in results:
                    fh.write("\t".join([str(i) for i in result]) + "\n")

    @staticmethod
    def _iter_from_file(f, batchsize):
        """ Yields groups of at most batchsize from a file. Results are lists of tuples """
        with open(f) as IN:
            lines_in_batch = 0
            batch = []
            for line in IN:
                batch.append(line[:-1].split("\t"))
                
                lines_in_batch += 1

                if lines_in_batch == batchsize:
                    yield batch
                    
                    batch = []
                    lines_in_batch = 0
                

    @staticmethod
    def _iter_sql_results(cursor, batchsize):
        """ Iterates over SQL results yielding batches of size=batchsize """
        while True:
            results = cursor.fetchmany(batchsize)
            if not results:
                break
            else:
                yield results

    def find_conserved_kmers(self, genomes=None, fuzzy=False):
        """ Finds conserved kmers in all (default) or a specific set of genomes (specified by genomes argument) """
        
        # make a temporary table with only the genomes we want to include
        # If genomes not specified, this table will be an exact clone on Genomes
        conn.execute("""CREATE TEMPORARY TABLE SubGen (
                        id INTEGER PRIMARY KEY NOT NULL,
                        name STRING
                        );
                        """)
        if genomes:
            conn.execute("""INSERT INTO SubGen (id, name)
                            SELECT id, name FROM Genomes
                            WHERE Genomes.id IN ({ph})
                            """.format(ph=",".join(["?"]*len(genomes))), genomes)
        else:
            conn.execute("""INSERT INTO SubGen (id, name)
                            SELECT id, name FROM Genomes
                            """)
 


        if fuzzy:

            conn.execute("SELECT * FROM FuzzyKmers WHERE (SELECT COUNT(DISTINCT Contigs.genome) FROM LOCATIONS JOIN Contigs ON Contigs.id = Locations.contig JOIN KmerToFuzzy ON KmerToFuzzy.fuzzy = FuzzyKmers.id WHERE Locations.kmer = KmerToFuzzy.kmer) = 2")

        else:

            """
            Here, we want to check, for each kmer, if that kmer is present in all genomes.

            SQL statement broken down:

            # select everything from the Kmers records that match
            SELECT * FROM Kmers 
                WHERE 

                # count unique genomes in the set of contigs generated by the following statements
                (SELECT COUNT(DISTINCT Contigs.genome)
            
                # select from the locations table
                FROM Locations
            
                # join with contigs to have access to genome data per location
                JOIN Contigs
                    ON Contigs.id=Locations.contig
            
                # right join with subgen to limit locations to entires within the genomes we want to consider
                RIGHT JOIN SubGen
                    ON SubGen.id = Contigs.genome

                # Link all the locations to their entry in Kmers
                WHERE Locations.kmer=Kmers.id)

            # check if num unique genomes == total genomes in the subset
            = SELECT COUNT(*) FROM SubGen;
            """
            conn.execute("""
                    SELECT * FROM Kmers 
                    WHERE 
                        (SELECT COUNT(DISTINCT Contigs.genome) 
                            FROM Locations 
                            JOIN Contigs 
                                ON Contigs.id=Locations.contig
                            RIGHT JOIN SubGen 
                                ON SubGen.id = Contigs.genome
                            WHERE Locations.kmer=Kmers.id) 
                        = SELECT COUNT(*) FROM SubGen;
                    """).fetchall() 


class KCounter(multiprocessing.Process):

    COMPLEMENT_MAP = {"A": "T", "T": "A", "G": "C", "C": "G"}

    def __init__(self, name, k, input_queue, results_queue):
        super().__init__()
        self.name = name
        self.k = k
        self.input_queue = input_queue
        self.results_queue = results_queue

        self.kmers = {}
        self.locations = []

    @classmethod
    def _reverse_complement(cls, kmer):
        return ''.join([cls.COMPLEMENT_MAP[nt] for nt in reversed(kmer)])

    @staticmethod
    def _kmer_to_id(kmer):
        """ Returns a base10 id from a base4 representation coding the nucleotides """
        # statements are chained together like this because this timed as the fastest way
        id_str = kmer.replace("A", "0").replace("C", "1").replace("G", "2").replace("T", "3")
        return int(id_str, 4)


    def run(self):
        LOG.debug("Starting {}".format(self.name)) 
        # wait for kmers to count
        while True:

            itm = self.input_queue.get()

            if type(itm) == tuple:
                genome, contig, start, seq = itm
                self.count_kmers(genome, contig, start, seq)
                self.dump_kmers_and_locations()
            elif itm == "STOP":
                LOG.debug("Returning from {}".format(self.name))
                return

            else:
                print("Not tuple!")



    def count_kmers(self, genome, contig, start, seq):
        for indx in range(len(seq)+1 - self.k):
            #kmers_counted += 1
            kmer = seq[indx:indx+self.k]

            # check if there is an N in the kmer and discard it if so
            if "N" in kmer:
                continue
               
            # choose the < of the kmer and its RC
            canonical = sorted((kmer, self._reverse_complement(kmer)))[0]
        
            # set strand based on whether or not rc was used
            if canonical == kmer:
                strand = "1"
            else:
                strand = "-1"

            kmer_id = self._kmer_to_id(kmer)
            # add the location and the kmer 
            self.locations.append((contig, start+indx, strand, kmer_id))
            self.kmers[kmer_id] = kmer
        return

    def dump_kmers_and_locations(self):
        """ Dumps the kmers and locations to the output queue and clears the variables from memory """
        # convert to kmer dict to a list of tuples
        kmer_dump = [(k, v) for k, v in self.kmers.items()]

        self.results_queue.put((kmer_dump, self.locations))

        # reset the variables
        self.locations = []
        self.kmers = {}
        

#=====================================================

class KmerCounter(object):

    # set max mem to 38G in mbytes
    MAX_MEM = 60 * 2 ** 10

    COMPLEMENT_MAP = {"A": "T", "T": "A", "G": "C", "C": "G"}

    def __init__(self, mismatches=1, stable_start=2, stable_end=2):
        self.mismatches = mismatches

        # stables are # bases at beginning and end that cannot be variable
        self.stable_start = stable_start
        self.stable_end = stable_end

        # make partitions based on stable_start
        self.queues = {}
        self.partitions = {}
        self.results_queue = multiprocessing.Queue()
        for kmer in self._recurse_all_kmers(self.stable_start):
            self.queues[kmer] = multiprocessing.Queue()
            self.partitions[kmer] = KmerPartition(kmer, self.queues[kmer], self.results_queue)
            self.partitions[kmer].start()

        # set genome id to starting value
        self.genome_id = 0

    @classmethod
    def _recurse_all_kmers(cls, k, curseq="", collected_seqs=[]):
        """ Generates all kmers to avoid having to store a large array """
        if len(curseq) == k:
            yield curseq
        else:
            for nt in ["A", "C", "G", "T"]:
                for kmer in cls._recurse_all_kmers(k, curseq + nt, collected_seqs):
                    yield kmer
 
    @classmethod
    def _reverse_complement(cls, kmer):
        return ''.join([cls.COMPLEMENT_MAP[nt] for nt in reversed(kmer)])


    def count_kmers(self, fasta_f, k):
        LOG.info("Counting {}-mers for {}".format(k, fasta_f))

        fasta_name = os.path.splitext(os.path.basename(fasta_f))[0]
        if os.path.isfile(fasta_name + ".counts.txt"):
            LOG.info("Found existing counts file for {}".format(fasta_name))
            kmer_counts = fasta_name + ".counts.txt"

        else:
            # count kmers using jellyfish
            jf_indx = jellyfish.count_kmers(fasta_f, k, prefix=fasta_name)
            kmer_counts = jellyfish.dump_kmer_counts(jf_indx, output_f=fasta_name + ".counts.txt")

        # read in the counts file
        num_processed = 0
        with open(kmer_counts, 'r') as IN:
            for line in IN:
                kmer, count = line.rstrip().split("\t")

                # send the kmer to the correct partition for processing
                self.partitions[kmer[:self.stable_start]].add_kmer(kmer, count)
                
                
                num_processed += 1
                #if num_processed == 100000:
                    #sys.exit()

        [print((par.name, len(par.kmers))) for par in self.partitions.values()]

        conserved_kmers = []
        for par in self.partitions.values():
            print((par.name, len(par.kmers)))
            par.update_kmer_file()
            conserved_kmers += par.get_conserved_kmers()

    def count_kmers_custom_single(self, fasta_f, k):
    
        self.genome_id += 1
        genome_str = str(self.genome_id)

        kmers_counted = 0

        contig_id = 0
        for seq in self._parse_fasta(fasta_f):
            contig_id += 1
            gen_ctg_str = ";".join((genome_str, str(contig_id)))

            # count kmers in a sliding window
            # must add 1 to the length because range has an exclusive end
            for indx in range(len(seq)+1 - k):
                kmers_counted += 1
                kmer = seq[indx:indx+k]

                # check if there is an N in the kmer and discard it if so
                if "N" in kmer:
                    continue
                
                self.partitions[kmer[:self.stable_start]].add_kmer(kmer, ";".join((gen_ctg_str, str(indx))))

                if kmers_counted % 10000 == 0:
                    print(kmers_counted)

        for partition in self.partitions.values():
            partition.update_kmer_file()



    def count_kmers_custom(self, fasta_f, k):
        """ Count kmers and their reverse complements """
    
        self.genome_id += 1
        genome_str = str(self.genome_id)

        kmers_counted = 0

        contig_id = 0
        for seq in self._parse_fasta(fasta_f):
            contig_id += 1
            contig_str = str(contig_id)
            #gen_ctg_str = ";".join((genome_str, str(contig_id)))

            # count kmers in a sliding window
            # must add 1 to the length because range has an exclusive end
            for indx in range(len(seq)+1 - k):
                kmers_counted += 1
                kmer = seq[indx:indx+k]

                # check if there is an N in the kmer and discard it if so
                if "N" in kmer:
                    continue
               

                # choose the < of the kmer and its RC
                canonical = sorted((kmer, self._reverse_complement(kmer)))[0]
        
                # check if the rc was used and reindex if not
                if canonical == kmer:
                    strand = "1"
            
                else:
                    strand = "-1"
 

                self.partitions[canonical[:self.stable_start]].kmer_queue.put((canonical, ";".join((genome_str, contig_str, str(indx), strand))))

                if kmers_counted % 10000 == 0:
                    # check if the queues are getting too big, or memory usage too high
                    for p in self.partitions.values():
                        if p.kmer_queue.qsize() > 1000:
                            time.sleep(2)

                    # check memory usage
                    main = psutil.Process(os.getpid())

                    # get mem of main process in MB
                    total_mem = main.get_memory_info()[0] / 2 ** 20
                    for partition in self.partitions.values():
                        ps = psutil.Process(partition.pid)
                        total_mem += ps.get_memory_info()[0] / 2 ** 20

                    # do an emergency dump if we are at 60% of max memory
                    # 60% because the writing to file step can add an extra 4% current mem per proc (16)
                    # and we don't want to exceed max
                    if total_mem >= .6 * self.MAX_MEM:
                        LOG.info("Memory usage = {}G... doing emergency dump.".format(total_mem / 2 ** 10))
                        for p in self.partitions.values():
                            p.kmer_queue.put("DUMP")
                    #print(kmers_counted)

        for p in self.partitions.values():
            p.kmer_queue.put("DUMP")


    def get_conserved_kmers_single(self, num_genomes):
        conserved_kmers = []
        for p in self.partitions.values():
            conserved_kmers += p.get_conserved_kmers()

        return conserved_kmers

    def get_conserved_kmers(self, num_genomes):
        for p in self.partitions.values():
            p.kmer_queue.put("CONSERVED{}".format(num_genomes))

        active_jobs = [p for p in self.partitions.values()]
        while active_jobs:
            try:
                yield(self.results_queue.get(True, 5))
            except queue.Empty:
                new_active_jobs = []
                for job in active_jobs:
                    if job.is_alive():
                        new_active_jobs.append(job)
                    else:
                        #job.kmer_queue.join_thread()
                        job.join()
                        job.kmer_queue.close()

                active_jobs = new_active_jobs
                time.sleep(5)

        LOG.debug("All threads joined!")

        self.results_queue.close()


    def _parse_fasta(self, fasta_f):
        """ Yields sequences from a fasta file """

        with open(fasta_f, 'r') as IN:
            for record in SeqIO.parse(IN, 'fasta'):
                yield str(record.seq)


class KmerPartition(multiprocessing.Process):
    """ A free-standing class (for multithreading purposes) that processes a subset of kmers 
    
    Responsible for holding kmer information and dumping it to a file when it gets the signal.

    I think the best way to hold kmer information is in a trie. 
    """

    # binary here is merely idealistic because the binary gets converted to ints
    #encode_map = {"A": 0b00 "a": 0b00, "C": 0b01, "c": 0b01, "G": 0b10, "g": 0b10, "T": 0b11, "t":0b11}
    encode_map = {"A": "00", "C": "01", "G": "10", "T": "11"}


    def __init__(self, name, kmer_queue, results_queue):
        super().__init__()
        self.name = name
        self.kmer_queue = kmer_queue
        self.results_queue = results_queue

        # set up the storage file
        self.out_dir = "."
        self.file = "/".join([self.out_dir, self.name + "_data_dump.txt"])
        
        # set up some vars to use later
        self.kmers = {}
    
        self.mismatches = 1
        self.stable_start = 2
        self.stable_end = 2
        self.k = 20


    def run(self):
        LOG.debug("Starting thread {}".format(self.name)) 
        # wait for kmers to count
        while True:
            itm = self.kmer_queue.get()

            if type(itm) is tuple:
                self.add_kmer(kmer=itm[0], position=itm[1])

            else:

                if itm == "DUMP":
                    self.update_kmer_file()

                elif itm.startswith("CONSERVED"):
                    # split the number of genomes from the itm
                    self.get_conserved_kmers(int(itm.split("D")[1]))
                    return

        return

    @classmethod
    def _recurse_all_groups(cls, k, curseq="", mismatches=0):
        """ Generates all kmer groups for a given length k """
        if len(curseq) == k:
            yield curseq

        # don't allow N's in stable regions
        elif len(curseq) <= self.stable_start:
            pass

            for nt in ["N", "A", "C", "G", "T"]:
                for kmer in cls._recurse_all_kmers(k, curseq + nt):
                    yield kmer
 

    def _digest_kmer(self, kmer, current_seq="", mismatches=0):
        """ This is a slow point of my program right now. I need to speed up. """
        # check if we are to the end
        if len(current_seq) >= self.k - self.stable_end:
            # not sure why I can't just return rather than yield, then return
            # might have something to do with the fact that some will be yielded and others will be returned
            yield "".join((current_seq, kmer))
            return
        
        # skip if we are still in the stable start
        elif len(current_seq) < self.stable_start:
            pass

        else:
            
            # check if the current base can be fuzzy (not too many previous errors)
            if mismatches < self.mismatches:
                # run the ambiguous result
                for result in self._digest_kmer(kmer[1:], current_seq+"N", mismatches+1):
                    yield result
                    break
 
            # base cannot be fuzzy, already enough errors
            else:
                # once again, not sure why I can't just return this...
                yield "".join((current_seq, kmer))
                return
        
        # run the next iteration if not returned
        for result in self._digest_kmer(kmer[1:], current_seq+kmer[0], mismatches):
            #print(("result", result))
            yield result

    def add_kmer(self, kmer, position):

        for group in self._digest_kmer(kmer):
            try:
                self.kmers[group].append(position)
            except KeyError:
                self.kmers[group] = [position]


    ###
    # Checking analyzing stored files
    ###
    def get_conserved_kmers(self, num_genomes):
        """ Look for kmers that have a position from each genomes. """
       
        LOG.debug("Finding conserved kmers for {}...".format(self.name))
        conserved_kmers = []
        genomes_not_found = {}
        with open(self.file, 'r') as IN:
            for line in IN:
                elems = line[:-1].split("\t")
                kmer = elems[0]
                locations = elems[1:]
    
                # let's assume that most kmers will NOT be found in all genomes
                # so, I don't want to make objects (slow) until I know they are conserved. 
                # this will be slower for kmers that ARE in all genomes
    
                # first check if there are enough locations -- this will remove most kmers
                if len(locations) >= num_genomes:
                    
                    # now check that each genome is represented and build a location list
                    genomes_found = [False]*num_genomes
                    for location in locations:

                        try:
                            genome = int(location.split(";")[0]) - 1
                        except:
                            print(locations)
                            raise
                        try:
                            genomes_found[genome].append(location)
                        except AttributeError:
                            genomes_found[genome] = [location]
                        
                    not_found = []
                    for indx, pos in enumerate(genomes_found):
                        if pos is False:
                            not_found.append(indx)
                    
      
                    if not_found:
                        # build a histogram like data structure with kmers not found
                        try:
                            genomes_not_found[len(not_found)].append(kmer)
                        except KeyError:
                            genomes_not_found[len(not_found)] = [kmer]
                    else:
                        # if all genomes found, convert location into object and add to list of conserved regions
                        loc_objs = []
                        for location in locations:
                            genome, contig, start, strand = location.split(";")
                            loc_objs.append(Location(genome, contig, start, strand))
                        conserved_kmers.append((kmer, loc_objs))

        LOG.debug("Found {} conserved kmers - {}".format(len(conserved_kmers), self.name))
        for kmer in conserved_kmers:
            self.results_queue.put(kmer)
        return conserved_kmers


    def expand_and_update_kmer_file(self):
        """ Loads/creates a kmer file and writes all the current kmer information to it. """

        LOG.debug("Updating kmer file for {}".format(self.file))

        tmp_file = self.file + "_tmp"

        # iterate through the input file and stored kmers and dump ordered results
        with open(self.file, 'r') as IN, open(tmp_file, 'w') as OUT:

            stored_line = IN.readline()[:-1]
           

            for kmer in sorted(self.kmers):
                # make sure something is on the line; assume file over if not
                if stored_line:
                    #stored_kmer = int(stored_line.split("\t")[0])
                    stored_kmer = stored_line.split("\t")[0]
                else:
                    line = "\t".join([str(kmer), "\t".join(self.kmers[kmer])])
                    OUT.write(line + "\n")
                    continue

                
                if kmer < stored_kmer:      # stored_kmer comes after this one, add new line
                    # write kmer and all locations
                    line = "\t".join([str(kmer),"\t".join(self.kmers[kmer])])

                elif kmer == stored_kmer:   # stored_kmer is this one, append to line
                    line = "\t".join([stored_line, "\t".join(self.kmers[kmer])])
                    stored_line = IN.readline()[:-1]

                else:   # stored_kmer is less than this one; need to read some more lines
                    # read in some more lines until the line is greater than 
                    line = stored_line
                    while kmer > stored_kmer:
                        stored_line = IN.readline()[:-1]

                        # check for the end of the file and break if found
                        if not stored_line:
                            break

                        OUT.write(line + "\n")
                        stored_kmer = stored_line.split("\t")[0]

                        if kmer == stored_kmer:
                            line = "\t".join([stored_line, "\t".join(self.kmers[kmer])])
                            stored_line = IN.readline()[:-1]
                            # while loop will break at the end of this
                        elif kmer < stored_kmer:    # we have overshot, add the kmer from memory
                            line = "\t".join([str(kmer), "\t".join(self.kmers[kmer])])
                            # while loop will break at the end of this

                        else:   # there is still more looping over file needed to be done
                            line = stored_line


                OUT.write(line + "\n")

        os.rename(tmp_file, self.file)
        self.kmers = {}


    def merge_kmers_with_file(self, fh):
        """ Yields lines to write to an output file """
        kmer_itr = iter(sorted(self.kmers))
        file_itr = iter(fh.readline, "")

        kmer = next(kmer_itr)
        line = next(file_itr)[:-1]
    
        while True:
            if kmer < line.split("\t")[0]:
                #print("<")
                yield "\t".join([str(kmer),"\t".join(self.kmers[kmer])])

                try:
                    kmer = next(kmer_itr)
                except StopIteration:
                    # no more amps, yield the rest of the file
                    yield line
                    while True:
                        yield next(file_itr)[:-1]

            elif kmer == line.split("\t")[0]:
                #print("=")
                yield "\t".join([line, "\t".join(self.kmers[kmer])])

                try:
                    kmer = next(kmer_itr)
                except StopIteration:
                    # no more kmers, yield the rest of the file
                    while True:
                        yield next(file_itr)[:-1]

                try:
                    line = next(file_itr)[:-1]
                except StopIteration:
                    while True:
                        kmer = next(kmer_itr)
                        yield "\t".join([str(kmer),"\t".join(self.kmers[kmer])])

            else:
                #print(">")
                yield line

                try:
                    line = next(file_itr)[:-1]
                except StopIteration:
                    yield "\t".join([str(kmer),"\t".join(self.kmers[kmer])])

                    while True:
                        kmer = next(kmer_itr)
                        yield "\t".join([str(kmer),"\t".join(self.kmers[kmer])])

    def update_kmer_file(self):

        if os.path.isfile(self.file):
            tmp_f = self.file + "_tmp"
            with open(self.file, 'r') as IN, open(tmp_f, 'w') as OUT:
                for line in self.merge_kmers_with_file(IN):
                    #print((line,))
                    #print()
                    OUT.write(line + "\n")

            os.rename(tmp_f, self.file)
        else:
            # file doesn't already exist just write all kmers
            with open(self.file, 'w') as OUT:
                for kmer in sorted(self.kmers):
                    OUT.write("\t".join([str(kmer),"\t".join(self.kmers[kmer])]) + "\n")

        self.kmers = {}


class Amplicon(object):

    amplicon_index = 0

    def __init__(self, name):
        self.name = name
        self.locations = []

        Amplicon.amplicon_index += 1
        self.index = Amplicon.amplicon_index

    def __str__(self):
        return self.name

    def __cmp__(self, other):
        return cmp(self.name, other.name)

    def __lt__(self, other):
        return self.name < other.name

    def __eq__(self, other):
        return self.name == other.name

    def __gt__(self, other):
        return self.name > other.name

    def _get_stored_repr(self, base=None):
        """ Generates a representation of the object that can be loaded later, optionally just adds locations to a base string """

        if base is None:
            #print(self.locations)
            return "\t".join([self.name] + [":".join((str(loc1), str(loc2))) for loc1, loc2 in self.locations])
        else:
            #print(self.locations)
            #print(base)
            return "\t".join([base] + [":".join((str(loc1), str(loc2))) for loc1, loc2 in self.locations])

    @classmethod
    def from_stored_repr(cls, rep):
        elems = rep.split("\t")
        amp = cls(elems[0])

        for loc in elems[1:]:
            l1, l2 = loc.split(":")

            amp.add_location_pair(Location.from_str(l1), Location.from_str(l2))

        return amp

    @classmethod
    def merge_amplicons_with_file(cls, amps, fh):
        """ Yields lines to write to an output file """
        amp_itr = iter(amps)
        file_itr = iter(fh.readline, "")

        amp = next(amp_itr)
        line = next(file_itr)[:-1]
    
        while True:
            #print((amp.name, line.split("\t")[0]))
            if amp.name < line.split("\t")[0]:
                #print("<")
                yield amp._get_stored_repr()

                try:
                    amp = next(amp_itr)
                except StopIteration:
                    # no more amps, yield the rest of the file
                    yield line
                    while True:
                        yield next(file_itr)[:-1]

            elif amp.name == line.split("\t")[0]:
                #print("=")
                yield amp._get_stored_repr(line)
                try:
                    amp = next(amp_itr)
                except StopIteration:
                    # no more amps, yield the rest of the file
                    while True:
                        yield next(file_itr)[:-1]

                try:
                    line = next(file_itr)[:-1]
                except StopIteration:
                    while True:
                        yield next(amp_itr)._get_stored_repr()

            else:
                #print(">")
                yield line

                try:
                    line = next(file_itr)[:-1]
                except StopIteration:
                    yield amp._get_stored_repr()
                    while True:
                        yield next(amp_itr)._get_stored_repr()

    @classmethod
    def update_amplicon_file(cls, amplicons, f):

        if os.path.isfile(f):
            tmp_f = f + "_tmp"
            with open(f, 'r') as IN, open(tmp_f, 'w') as OUT:
                for line in cls.merge_amplicons_with_file(sorted(amplicons), IN):
                    #print(line)
                    #print()
                    OUT.write(line + "\n")

            os.rename(tmp_f, f)
        else:
            # file doesn't already exist just write all amplicons
            with open(f, 'w') as OUT:
                for amp in sorted(amplicons):
                    OUT.write(amp._get_stored_repr() + "\n")

            
    def update(self, other):
        """ Updates an amplicon with more locations """
        self.locations = self.locations + other.locations

    def add_location_pair(self, loc1, loc2):
       self.locations.append((loc1, loc2))


class Location(object):
    """
    A location of a kmer indexed by genome, contig, start, and strand 
    
    Uses slots to be a lightweight struct.

    Fun fact, this is actually smaller than a tuple with the same info (this takes 72 bytes, tup is 80) according to sys.getsizeof()
    """
    #__slots__ = ["genome", "contig", "start", "strand"]

    def __init__(self, genome, contig, start, strand):
        self.genome = int(genome)
        self.contig = int(contig)
        self.start = int(start)
        self.strand = int(strand)

    @classmethod
    def from_str(cls, string):
        """ Returns a location object from a ';' sparated string representation """
        genome, contig, start, strand = string.split(";")
        return Location(genome, contig, start, strand)
    
    def __str__(self):
        return ";".join((str(self.genome), str(self.contig), str(self.start), str(self.strand)))

#
## Find kmers that are the right distance apart
#



def get_properly_spaced_pairs(conserved_kmers, k=20, max_dist=400, min_dist=100):

    # build partitions by contig because we know that to be the right distance, it must be on the same contig
    # this just reformats the data for a quicker searching algorithm
    ppp = {}
    for kmer, locations in conserved_kmers:
        for location in locations:
            
            try:
                ppp[(location.genome, location.contig)].append((kmer, location))
            except KeyError:
                ppp[(location.genome, location.contig)] = [(kmer, location)]

    # this part could be split up over multiple processes
    LOG.info("Processing {} conserved kmer partitions...".format(len(ppp)))
    results_queue = multiprocessing.Queue()
    bins_processed = 1
    for data_bin in ppp.values():
        LOG.debug("Processing bin {}...".format(bins_processed))
        bins_processed += 1


        for itm in _find_properly_spaced(results_queue, data_bin, k, max_dist, min_dist):
            yield itm
"""
    while True:
        if results_queue.qsize() > 0:
            yield results_queue.get()

        else:
            break
"""

class AmpliconFinder(object):
    """ Organizes all the code associated with going from a list of conserved kmers to a finished product of amplicons """

    def __init__(self, threads=7, k=20, max_dist=400, min_dist=100, max_mem=60):
        self.threads = threads
        self.k = k
        self.max_dist = max_dist
        self.min_dist = min_dist

        # attributes to be set later
        self.worker_pool = None


    def create_worker_pool(self):
        """ Makes a worker pool before loading all conserved kmers into memory so the conserved kmers are not duplicated. """
        self.worker_pool = multiprocessing.Pool(self.threads, self.worker_get_properly_spaced, (self,)) 

        # we can go ahead and close it because we are communicating with queues not tasks
        self.worker_pool.close()

        
    def find_kmer_pairs(conserved_kmers):
        """ 
        Writes a file of potential amplicons and returns the path to it. 
    
        Decides whether it would be best to look for amplicons that include all genomes or amplicons that include some minimum number of genomes.
        """

        if len(conserved_kmers) > 500:
            # want to do the limited search
            pass
        else:
            # want to do the full search
            pass

    def _find_kmer_pairs_in_all(conserved_kmers):
        """ Branch of find_kmer_pairs that looks for kmers that are present in all genomes """
        pass

    def _find_kmer_pairs_in_genomes(conserved_kmers, min_number_genomes):
        """ 
        Branch of find_kmer_pairs that looks for kmers that are present in some number of genomes.
        
        The Problem:

        Given a set of kmers with at least one location in each genome, find all possible kmer pairs that are spaced the proper distance apart in at least N genomes. 

        ===============
        Naiive Solution:

        All-by-all pair each kmer with each other kmer and count how many genomes have properly spaced locations. Compare this count to N.

        This will not work because sometimes genome sets have several million conserved kmers (resulting in more than 1 trillion pairwise comparisons.

        =============
        Optimizations:

        1. Can take advantage of small proper distance window by sorting kmers and comparing each kmer to those kmers located within a proper window. This allows me to compare each kmer to at most 300 (400 max range, 100 min range) * all possible N variants

        2. Can stop looking when # properly spaced genomes + # remaining genomes < # genomes required

        ========
        Strategy:

        Current:
        1. Split kmers into locations per contigs
        2. Send contigs to threads and threads return good hits
        3. Make good hits into amplicons and write to file when mem usage too high


        Alternate:
        1. Split kmers into their locations per genome (and per contig just to save a step later)
        2. Give dict of locations corresponding to a genome to each thread
        3. Each thread calculates all properly spaced pairs and writes to a file (some sort of memory monitoring should occur, either by main or each thread should be assigned memory)
        4. Multi-file merge sort into a single file

        This schema has the advantage that intermittant writting to disk is done threadded and also main is not a big relay system that all threads must pass through. Also, each result need not be pickled to be put in a queue.

        This schema will be much faster iff the main thread is the primary limiting factor. (I suspect it is) 

        ---- Neither of these take advantage of Optimization #2.


        To use optimization 2, we have to be kmer focused rather than genome focused because we need to be able to know how many genomes have that kmer pair in a proper spacing.

        Op2 Solution:
        1. Main thread iterates over all kmer pairs and puts them in a queue.
        2. Workers put only ones that pass in an output queue

        This takes no advantage of Op1 -- which may end up being the one that saves more comparisons.

        """
        pass


def create_worker_pool(threads, input_queue, output_queue, k, max_dist, min_dist):
    """ Makes a worker pool before loading all conserved kmers into memory so the conserved kmers are not duplicated. """
    pool = multiprocessing.Pool(threads, worker_get_properly_spaced, (input_queue, output_queue, k, max_dist, min_dist)) 
    pool.close()

    return pool


def pooled_find_kmer_pairs(pool, conserved_kmers, input_queue, output_queue, threads, max_mem=60):
    """ Even after dumping, memory usage doesn't go down as much as it should """

    # build partitions by contig because we know that to be the right distance, it must be on the same contig
    # this just reformats the data for a quicker searching algorithm
    #objgraph.show_growth()
    ppp = {}
    for kmer, locations in conserved_kmers:
        for location in locations:
            
            try:
                ppp[(location.genome, location.contig)].append((kmer, location))
            except KeyError:
                ppp[(location.genome, location.contig)] = [(kmer, location)]

    # add a start tag to the queue
    for _ in range(threads):
        input_queue.put("START")

    # add all bins to the queue
    for data_bin in ppp.values():
        input_queue.put(data_bin)


    LOG.info("Processing {} conserved kmer partitions...".format(len(ppp)))
    #tr.print_diff()

    # get all the results
    results_processed = 0
    amplicons = {}
    num_finished = 0
    while True:

        # check memory every so often and dump if necessary
        if results_processed != 0 and results_processed % 500000 == 0:
            #LOG.debug("{}G memory used.".format(get_memory_usage()))
            if get_memory_usage() >= .7 * max_mem:
                LOG.debug("{}G memory used ({} results processed) -- Dumping amplicons to file.".format(get_memory_usage(True), results_processed)) 
                Amplicon.update_amplicon_file(list(amplicons.values()), "amplicons_data_dump.txt")
                amplicons = {}

                #objgraph.show_growth()
                #tr.print_diff()

                LOG.debug("{}G memory used after dumping.".format(get_memory_usage()))

            LOG.info("{} partitions remaining".format(input_queue.qsize()))

        # get a result from the queue waiting up to 5 seconds
        try:
            result = output_queue.get(True, 5)
        except queue.Empty:
            if num_finished == threads:
                LOG.debug("All threads done. Joining.")
                # Try to join the threads
                pool.join()
                LOG.debug("{} Results (locations) processed.".format(results_processed))
                break
            else:
                continue
                
            # need some sort of check if the work is done before exiting

        if result == "FINISHED":
            num_finished += 1
            continue

        key, loc1, loc2 = result

        try:
            amplicons[key].add_location_pair(loc1, loc2)
        except KeyError:
            amp = Amplicon(key)
            amp.add_location_pair(loc1, loc2)
            amplicons[key] = amp
        
        results_processed += 1

    # final pickle
    Amplicon.update_amplicon_file(list(amplicons.values()), "amplicons_data_dump.txt")

    return "amplicons_data_dump.txt"

def pooled_get_properly_spaced_pairs(conserved_kmers, k=20, max_dist=400, min_dist=100, max_mem=60):

    # build partitions by contig because we know that to be the right distance, it must be on the same contig
    # this just reformats the data for a quicker searching algorithm
    #objgraph.show_growth()
    ppp = {}
    for kmer, locations in conserved_kmers:
        for location in locations:
            
            try:
                ppp[(location.genome, location.contig)].append((kmer, location))
            except KeyError:
                ppp[(location.genome, location.contig)] = [(kmer, location)]

    
    # make queue and add all entries to it
    input_queue = multiprocessing.Queue()
    for data_bin in ppp.values():
        input_queue.put(data_bin)


    LOG.debug("Calculating memory.")
    objs = muppy.get_objects()
    sum1 = summary.summarize(objs)
    summary.print_(sum1)

    LOG.debug("Done calculating memory.")

    print(dir())


    sys.exit()


    output_queue = multiprocessing.Queue(10000)
    threads = 7
    pool = multiprocessing.Pool(threads, worker_get_properly_spaced, (input_queue, output_queue, k, max_dist, min_dist)) 
    pool.close()

    LOG.info("Processing {} conserved kmer partitions...".format(len(ppp)))
    #tr.print_diff()

    # get all the results
    results_processed = 0
    amplicons = {}
    num_finished = 0
    while True:

        # check memory every so often and dump if necessary
        if results_processed != 0 and results_processed % 100000 == 0:
            #LOG.debug("{}G memory used.".format(get_memory_usage()))
            if get_memory_usage() >= .7 * max_mem:
                LOG.debug("{}G memory used ({} results processed) -- Dumping amplicons to file.".format(get_memory_usage(), results_processed)) 
                Amplicon.update_amplicon_file(list(amplicons.values()), "amplicons_data_dump.txt")
                amplicons = {}

                #objgraph.show_growth()
                #tr.print_diff()

                LOG.debug("{}G memory used after dumping.".format(get_memory_usage()))

            LOG.info("{} partitions remaining".format(input_queue.qsize()))

        # get a result from the queue waiting up to 5 seconds
        try:
            result = output_queue.get(True, 5)
        except queue.Empty:
            if num_finished == threads:
                LOG.debug("All threads done. Joining.")
                # Try to join the threads
                pool.join()
                LOG.debug("{} Results (locations) processed.".format(results_processed))
                break
            else:
                continue
                
            # need some sort of check if the work is done before exiting

        if result == "FINISHED":
            num_finished += 1
            continue

        key, loc1, loc2 = result

        try:
            amplicons[key].add_location_pair(loc1, loc2)
        except KeyError:
            amp = Amplicon(key)
            amp.add_location_pair(loc1, loc2)
            amplicons[key] = amp
        
        results_processed += 1

    # final pickle
    Amplicon.update_amplicon_file(list(amplicons.values()), "amplicons_data_dump.txt")

    return "amplicons_data_dump.txt"


def pooled_get_properly_spaced_pairs_option(conserved_kmers, k=20, max_dist=400, min_dist=100, max_mem=60):

    # build partitions by contig because we know that to be the right distance, it must be on the same contig
    # this just reformats the data for a quicker searching algorithm
    #objgraph.show_growth()
    ppp = {}
    for kmer, locations in conserved_kmers:
        for location in locations:
            
            try:
                ppp[(location.genome, location.contig)].append((kmer, location))
            except KeyError:
                ppp[(location.genome, location.contig)] = [(kmer, location)]

    
    # make queue and add all entries to it
    input_queue = multiprocessing.Queue()
    for data_bin in ppp.values():
        input_queue.put(data_bin)


    LOG.debug("Calculating memory.")
    objs = muppy.get_objects()
    sum1 = summary.summarize(objs)
    summary.print_(sum1)

    LOG.debug("Done calculating memory.")

    print(dir())


    sys.exit()


    output_queue = multiprocessing.Queue(10000)
    threads = 7
    pool = multiprocessing.Pool(threads, worker_get_properly_spaced, (input_queue, output_queue, k, max_dist, min_dist)) 
    pool.close()

    LOG.info("Processing {} conserved kmer partitions...".format(len(ppp)))
    #tr.print_diff()

    # get all the results
    results_processed = 0
    amplicons = {}
    num_finished = 0
    while True:

        # check memory every so often and dump if necessary
        if results_processed != 0 and results_processed % 100000 == 0:
            #LOG.debug("{}G memory used.".format(get_memory_usage()))
            if get_memory_usage() >= .7 * max_mem:
                LOG.debug("{}G memory used ({} results processed) -- Dumping amplicons to file.".format(get_memory_usage(), results_processed)) 
                Amplicon.update_amplicon_file(list(amplicons.values()), "amplicons_data_dump.txt")
                amplicons = {}

                #objgraph.show_growth()
                #tr.print_diff()

                LOG.debug("{}G memory used after dumping.".format(get_memory_usage()))

            LOG.info("{} partitions remaining".format(input_queue.qsize()))

        # get a result from the queue waiting up to 5 seconds
        try:
            result = output_queue.get(True, 5)
        except queue.Empty:
            if num_finished == threads:
                LOG.debug("All threads done. Joining.")
                # Try to join the threads
                pool.join()
                LOG.debug("{} Results (locations) processed.".format(results_processed))
                break
            else:
                continue
                
            # need some sort of check if the work is done before exiting

        if result == "FINISHED":
            num_finished += 1
            continue

        key, loc1, loc2 = result

        try:
            amplicons[key].add_location_pair(loc1, loc2)
        except KeyError:
            amp = Amplicon(key)
            amp.add_location_pair(loc1, loc2)
            amplicons[key] = amp
        
        results_processed += 1

    # final pickle
    Amplicon.update_amplicon_file(list(amplicons.values()), "amplicons_data_dump.txt")

    return "amplicons_data_dump.txt"


def worker_get_properly_spaced(input_queue, output_queue, k, max_dist, min_dist):

    start_signal = input_queue.get(True)
    if start_signal == "START":
        # short sleep to reduce the risk off this process getting another process' start signal
        time.sleep(3)
    else:
        raise ValueError("Data received before START signal. Possibly because of race condition.")

    while True:
        # get partitions from queue until it is empty then break
        try:
            data_bin = input_queue.get(True, 5)
        except queue.Empty:
            output_queue.put("FINISHED")
            break

        # sort kmers by position
        sorted_kmers = sorted(data_bin, key=lambda k: k[1].start)

        #print("# kmers: {}".format(len(sorted_kmers)))

        num_properly_spaced = 0
        for indx1, k1 in enumerate(sorted_kmers):
            # begin at the next kmer
            indx2 = indx1 + 1
            while indx2 < len(sorted_kmers):
                k2 = sorted_kmers[indx2]
                indx2 += 1

                # Find distance between start of k2 and the end of k1 (interprimer distance)
                # locations are starting values so need to add k to k1 to get the end of it
                # We know k2 is always >= k1 because we sorted it
                distance = k2[1].start - (k1[1].start + k) 

                # break at first one that fails max dist test
                if distance > max_dist:
                    break

                # skip if the distance is too small
                elif distance < min_dist:
                    continue

                # distance must be acceptable
                else:

                    # make sure they are on the same strand
                    if k1[1].strand == k2[1].strand:

                        key = "-".join(sorted((k1[0], k2[0])))
                        # puts a tuple of (key, loc1, loc2)
                        output_queue.put((key, k1[1], k2[1]))

def worker_get_properly_spaced_options(input_queue, output_queue, k, max_dist, min_dist):
    """ Worker thread that gets told if it should look for kmers amplicons present in ANY genome or only amplicons present in ALL genomes """

    method = input_queue.get(True)
    if method == "START_ALL" or method == "START_ANY":
        # short sleep to reduce the risk off this process getting another process' start signal
        time.sleep(3)
    else:
        raise ValueError("Data received before START signal. Possibly because of race condition.")

    while True:
        # get partitions from queue until it is empty then break
        try:
            data_bin = input_queue.get(True, 5)
        except queue.Empty:
            output_queue.put("FINISHED")
            break


        # Start looking for kmers that are present in ALL genomes
        # data_bin for this will be a kmer object 
        if method == "START_ALL":
            pass











        # sort kmers by position
        sorted_kmers = sorted(data_bin, key=lambda k: k[1].start)

        #print("# kmers: {}".format(len(sorted_kmers)))

        num_properly_spaced = 0
        for indx1, k1 in enumerate(sorted_kmers):
            # begin at the next kmer
            indx2 = indx1 + 1
            while indx2 < len(sorted_kmers):
                k2 = sorted_kmers[indx2]
                indx2 += 1

                # Find distance between start of k2 and the end of k1 (interprimer distance)
                # locations are starting values so need to add k to k1 to get the end of it
                # We know k2 is always >= k1 because we sorted it
                distance = k2[1].start - (k1[1].start + k) 

                # break at first one that fails max dist test
                if distance > max_dist:
                    break

                # skip if the distance is too small
                elif distance < min_dist:
                    continue

                # distance must be acceptable
                else:

                    # make sure they are on the same strand
                    if k1[1].strand == k2[1].strand:

                        key = "-".join(sorted((k1[0], k2[0])))
                        # puts a tuple of (key, loc1, loc2)
                        output_queue.put((key, k1[1], k2[1]))

def remove__find_properly_spaced(queue, data_bin, k=20, max_dist=400, min_dist=100):
    """ 
    Must look for the proper location.

    Must also be on the same strand???
    I think that's true but I should think about it more.
    """
    # sort kmers by position
    sorted_kmers = sorted(data_bin, key=lambda k: k[1].start)

    #print("# kmers: {}".format(len(sorted_kmers)))

    num_properly_spaced = 0
    for indx1, k1 in enumerate(sorted_kmers):
        # begin at the next kmer
        indx2 = indx1 + 1
        while indx2 < len(sorted_kmers):
            k2 = sorted_kmers[indx2]
            indx2 += 1

            # Find distance between start of k2 and the end of k1 (interprimer distance)
            # locations are starting values so need to add k to k1 to get the end of it
            # We know k2 is always >= k1 because we sorted it
            distance = k2[1].start - (k1[1].start + k) 

            # break at first one that fails max dist test
            if distance > max_dist:
                break

            # skip if the distance is too small
            elif distance < min_dist:
                continue

            # distance must be acceptable
            else:

                # make sure they are on the same strand
                if k1[1].strand == k2[1].strand:

                    # puts a tuple(first, second) of tuple(kmer, location)
                    #queue.put((k1, k2))
                    num_properly_spaced += 1
                    yield (k1, k2)

    LOG.debug("Found {} properly spaced pairs.".format(num_properly_spaced))
    

def write_amplicons_from_file(amplicon_f, fastas, k, out_dir, soft_limit=500, hard_limit=2000):
    """ Uses heuristics to write a set number of ampliconsto fastas """
    num_amplicons = 0
    num_locations = 0
    with open(amplicon_f, 'r') as IN:
        for line in IN:
            amp = Amplicon.from_stored_repr(line[:-1])

            for location in amp.locations:
                num_locations += 1

            num_amplicons += 1

    LOG.info("{} Total amplicons detected.".format(num_amplicons))
    LOG.info("{} Total locations detected.".format(num_locations))


    # shift to writing amplicons to file
    if os.path.isdir(out_dir):
        raise ValueError("{} already exists, cowardly refusing to overwrite.".format(out_dir))
    else:
        os.mkdir(out_dir)


    num_genomes = len(fastas)

    param_list = [
            (num_genomes, 0, 0),    # perfect matches only
            (num_genomes, 0),       # allow Ns in primer regions
            (num_genomes),          # just require a seq from each genome present
            (0),                    # require nothing
            ]

    # iterate over the parameters getting sets of amplicons
    amplicons = []
    hard_limit_omit = 0
    for params in param_list:
        for amp in subset_amplicon_file(amplicon_f, *params):

            # check if amp already in amplicons
            if amp in amplicons:
                continue
            else:
                if len(amplicons) >= hard_limit:
                    hard_limit_omit += 1
                else:
                    amplicons.append(amp)

        if hard_limit_omit:
            LOG.warning("Reached amplicon hard limit in the middle of a step. {} similar amplicons omitted.".format(hard_limit_omit))
            break

    LOG.info("Continuing processing with {} amplicons.".format(len(amplicons)))

    LOG.info("Making sequence shopping list...")
    shopping_list = {}
    for amplicon in amplicons:
        
        for location in amplicon.locations:
            first, second = location

            # there has GOT to be a better way to do this...
            try:
                shopping_list[first.genome][first.contig][(first, second)].append(amplicon)
            except KeyError:
                try:
                    shopping_list[first.genome][first.contig][(first, second)] = [amplicon]
                except KeyError:
                    try:
                        shopping_list[first.genome][first.contig] = {(first, second): [amplicon]}
                    except KeyError:
                        shopping_list[first.genome] = {first.contig: {(first, second): [amplicon]}}

    written_files = []
    for g_indx, fasta in enumerate(fastas, 1):
        if g_indx in shopping_list:
            with open(fasta) as IN:
                fasta_name = os.path.splitext(os.path.basename(fasta))[0]

                ### this part is getting down to a single region in the genome
                for c_indx, record in enumerate(SeqIO.parse(IN, "fasta"), 1):
                    if c_indx in shopping_list[g_indx]:
                        for location_pair in shopping_list[g_indx][c_indx]:
                            first, second = location_pair

                            # the current algorithm for finding amplicons ensures that position of first is always <= second, however, if it isn't the seq gotten from indexing will be blank (silently) so I want to make sure first is always before second
                            assert(first.start < second.start)
                            
                            # reverse complement if on rev strand
                            if first.strand == 1:
                                new_rec = record[first.start:second.start+k]
                            else:
                                new_rec = record[first.start:second.start+k].reverse_complement()

                            new_rec.id = "genome={f.genome};contig={f.contig};location=[{f.start}:{end}];strand={f.strand}".format(f=first, end=second.start + k)
                            new_rec.description = ""



                            ### this part is writing the region to each amplicon it is a part of
                            for amp in shopping_list[g_indx][c_indx][location_pair]:
                                out_path = "{}/{}_{}.fasta".format(out_dir, amp.index, amp.name)
                                written_files.append(out_path)
                                with open(out_path, 'a') as OUT:
                                    SeqIO.write(new_rec, OUT, "fasta")

    return written_files


def subset_amplicon_file(amplicon_f, min_genomes=0, max_duplicates=99, max_Ns=99):
    """ This could use to write unused amplicons to a file each time to narrow the search set with each iteration """
    with open(amplicon_f, 'r') as IN:
        for line in IN:
            elems = line[:-1].split("\t")

            #
            ## Tests I can do before making the amplicon
            #

            # quick check if amp has enough entries to cover each genome
            if len(elems) - 1 < min_genomes:
                continue
            
            # check the number of Ns
            if elems[0].count("N") > max_Ns:
                continue


            #
            ## Tests after making the amplicon
            #

            # rigorous num genomes test
            amplicon = Amplicon.from_stored_repr(line[:-1])
            genomes_represented = set()
            for location in amplicon.locations:
                genomes_represented.add(location[0].genome)

            if len(genomes_represented) < min_genomes:
                continue

            
            # num duplicates test
            num_duplicates = len(amplicon.locations) - len(genomes_represented)
            if num_duplicates > max_duplicates:
                continue

            yield amplicon


def remove__write_all_amplicons(amplicons, fastas, k, out_dir):
    """ 
    TODO: Need to do something with kmers that are duplicated in a genome.

    Maybe it would be best to take it into account in the stats section...

    Also need to figure out why RC isn't working correctly. All strands ate +1
    """

    if os.path.isdir(out_dir):
        raise ValueError("{} already exists, cowardly refusing to overwrite.".format(out_dir))
    else:
        os.mkdir(out_dir)

   
    #
    ## make a hierarchial list of amplicons to get to read each input FASTA once
    #
    LOG.info("Making sequence shopping list...")
    num_amplicons = 0
    shopping_list = {}
    for amplicon in amplicons.values():
        # make sure the amplicon represents more than one genome
        # this is a pre-check for adding to the shopping list
        genomes_represented = set()
        for location in amplicon.locations:
            genomes_represented.add(location[0].genome)

        if len(genomes_represented) == 1:
            continue

       
        
        num_amplicons += 1
        for location in amplicon.locations:
            first, second = location

            # there has GOT to be a better way to do this...
            try:
                shopping_list[first.genome][first.contig][(first, second)].append(amplicon)
            except KeyError:
                try:
                    shopping_list[first.genome][first.contig][(first, second)] = [amplicon]
                except KeyError:
                    try:
                        shopping_list[first.genome][first.contig] = {(first, second): [amplicon]}
                    except KeyError:
                        shopping_list[first.genome] = {first.contig: {(first, second): [amplicon]}}

    LOG.info("Writting {} amplicon files...".format(num_amplicons))
    written_files = []
    for g_indx, fasta in enumerate(fastas, 1):
        if g_indx in shopping_list:
            with open(fasta) as IN:
                fasta_name = os.path.splitext(os.path.basename(fasta))[0]

                ### this part is getting down to a single region in the genome
                for c_indx, record in enumerate(SeqIO.parse(IN, "fasta"), 1):
                    if c_indx in shopping_list[g_indx]:
                        for location_pair in shopping_list[g_indx][c_indx]:
                            first, second = location_pair

                            # the current algorithm for finding amplicons ensures that position of first is always <= second, however, if it isn't the seq gotten from indexing will be blank (silently) so I want to make sure first is always before second
                            assert(first.start < second.start)
                            
                            # reverse complement if on rev strand
                            if first.strand == 1:
                                new_rec = record[first.start:second.start+k]
                            else:
                                new_rec = record[first.start:second.start+k].reverse_complement()

                            new_rec.id = "genome={f.genome};contig={f.contig};location=[{f.start}:{end}];strand={f.strand}".format(f=first, end=second.start + k)
                            new_rec.description = ""



                            ### this part is writing the region to each amplicon it is a part of
                            for amp in shopping_list[g_indx][c_indx][location_pair]:
                                out_path = "{}/{}_{}.fasta".format(out_dir, amp.index, amp.name)
                                written_files.append(out_path)
                                with open(out_path, 'a') as OUT:
                                    SeqIO.write(new_rec, OUT, "fasta")

    return written_files

def get_amplicon_stats(amplicon_fastas, input_fastas, k):
    
    LOG.info("Calculating Amplicon stats...")

    # make a list of fasta names to be used in the matrix
    input_names = [os.path.splitext(os.path.basename(f))[0] for f in input_fastas]

    amplicon_stats = {}
    for fasta_f in amplicon_fastas:
        fasta_name = os.path.splitext(os.path.basename(fasta_f))[0]
        amplicon_stats[fasta_name] = {}

        # get the number of Ns from the name
        # fasta name will have the format $ampindx_$fwdprim-$revprim
        primers = fasta_name.split("_")[1]
        fwd_primer, rev_primer = primers.split("-")
       
        amplicon_stats[fasta_name]["fwd_N"] = str(fwd_primer.count("N"))
        amplicon_stats[fasta_name]["rev_N"] = str(rev_primer.count("N"))


        #
        ## check which seqs can be told apart
        #
        seqs = {}
        with open(fasta_f, 'r') as IN:
            for record in SeqIO.parse(IN, "fasta"):
                m = re.match(".*genome=(?P<genome>[^;]+).*strand=(?P<strand>[^;]+)", record.description)
                genome = m.group("genome")
                strand = m.group("strand")
                
                # store the amplified region (rc if necessary) 
                if strand == "1":
                    seq = record.seq[k:-k]
                else:
                    seq = record.seq[k:-k].reverse_complement()
                
                try:
                    seqs[genome].append(seq)
                except KeyError:
                    seqs[genome] = [seq]

        

        # Seqs can have multiple locations in a single amplicon -- flatten these
        flattened_seqs = {}
        for genome in seqs:
            for indx, seq in enumerate(seqs[genome]):
                flattened_seqs[(genome, indx)] = seq
        

        # Assign cluster for each seq
        cluster = 0
        cluster_counts = {}
        clusters = {}
        for seq1 in flattened_seqs:

            # skip seqs that already have cluster assigned
            if seq1 in clusters:
                continue
            else:
                cluster += 1
                clusters[seq1] = str(cluster)
                cluster_counts[cluster] = 1

            for seq2 in flattened_seqs:
                if seq1 is seq2:
                    continue

                # skip seqs that have already been processed
                elif seq2 in clusters:
                    continue

                else:
                    # check if the seqs are == 
                    if flattened_seqs[seq1] == flattened_seqs[seq2]:
                        clusters[seq2] = str(cluster)
                        cluster_counts[cluster] += 1

       
        # group clusters back by genomes
        duplicated_genomes = 0
        unique_genomes = set()     
        genome_clusters = {}
        for key in clusters:
            genome, index = key

            cluster = clusters[key]

            if cluster_counts[int(cluster)] == 1:
                unique_genomes.add(genome)
            try:
                genome_clusters[genome].append(clusters[key])
                duplicated_genomes += 1
            except KeyError:
                genome_clusters[genome] = [clusters[key]]
           
        amplicon_stats[fasta_name]["num_genomes"] = len(genome_clusters)
        amplicon_stats[fasta_name]["num_unique"] = len(unique_genomes)
        amplicon_stats[fasta_name]["num_duplicates"] = duplicated_genomes
        amplicon_stats[fasta_name]["genome_clusters"] = genome_clusters

    # sort the amplicons by num unique_genomes number (descending) and number of N's(ascending), and then num duplicates
    sorted_amplicons = sorted(amplicon_stats, key=lambda k: (-1 * amplicon_stats[k]["num_unique"], amplicon_stats[k]["fwd_N"] + amplicon_stats[k]["rev_N"], amplicon_stats[k]["num_duplicates"]))


    with open("amplicon_stats.txt", 'w') as STATS, open("amplicon_matrix.txt", 'w') as MATR:
        # write headers to both files
        STATS.write("\t".join(("amplicon", "num_genomes", "num_unique", "num_duplicates", "fwd_N", "rev_N")) + "\n")
        MATR.write("\t".join(["amplicon"] + input_names) + "\n")

        for amp in sorted_amplicons:
            write_stats = [amp] + [str(amplicon_stats[amp][k]) for k in ["num_genomes", "num_unique", "num_duplicates", "fwd_N", "rev_N"]]
            write_matr = [amp]
            for name in input_names:
                try:
                    write_matr.append(";".join(amplicon_stats[amp]["genome_clusters"][name]))
                except KeyError:
                    write_matr.append("0")

            STATS.write("\t".join(write_stats) + "\n")
            MATR.write("\t".join(write_matr) + "\n")


def get_memory_usage(verbose=False):
    """ Returns the memory usage of this process and all children """
    # check memory usage
    main = psutil.Process(os.getpid())

    total_mem = main.get_memory_info()[0] / 2 ** 30
    if verbose:
        LOG.debug("Main mem = {}".format(total_mem))
    for child in main.get_children():
        total_mem += child.get_memory_info()[0] / 2 ** 30
        if verbose:
            LOG.debug("    Thread usage {}".format(child.get_memory_info()[0] / 2 ** 30))

    return total_mem


def test_amplicons():
    k1 = ("AAAA-AAAA", Location(1, 1, 1, 1), Location(1, 1, 300, 1))
    k2 = ("AAAA-AAAA", Location(1, 2, 1, 1), Location(1, 5, 300, 1))
    k3 = ("CCCC-CCCC", Location(1, 3, 1, 1), Location(1, 2, 300, 1))
    k4 = ("AAAA-AAAA", Location(1, 4, 1, 1), Location(1, 1, 300, 1))
    k5 = ("GGGG-GGGG", Location(1, 5, 1, 1), Location(1, 1, 300, 1))

    amplicons = {}

    # add k1
    amp1 = Amplicon(k1[0])
    amp1.add_location_pair(k1[1], k1[2])
    amplicons[k1[0]] = amp1

    # add k2
    amplicons[k2[0]].add_location_pair(k2[1], k2[2])

    # add k3
    amp = Amplicon(k3[0])
    amp.add_location_pair(k3[1], k3[2])
    amplicons[k3[0]] = amp

    # dump to file
    Amplicon.update_amplicon_file(list(amplicons.values()), "test_file")
    amplicons = {}

    # add k4
    amp2 = Amplicon(k4[0])
    amp2.add_location_pair(k4[1], k4[2])
    amplicons[k4[0]] = amp2

    # add k5
    amp = Amplicon(k5[0])
    amp.add_location_pair(k5[1], k5[2])
    amplicons[k5[0]] = amp

    Amplicon.update_amplicon_file(list(amplicons.values()), "test_file")


    print("before update1")
    for pair in amp1.locations:
        print(pair[0].contig)

    print("before update2")
    for pair in amp2.locations:
        print(pair[0].contig)



    print("after update")
    amp1.update(amp2)
    for pair in amp1.locations:
        print(pair[0].contig)

    print("___________________________________")

    # read the test_file
    with open("test_file", 'rb') as IN:
        while True:
            amp = pickle.load(IN)
            print((amp.name, amp.locations, len(amp.locations)))
            for pair in amp.locations:
                print(pair[0].contig)



def main(fastas, k):
    #test_amplicons()
    kcounter = KmerCounter()
   
    # count kmers
    for fasta in fastas:
        LOG.info("Indexing k-mers in {}".format(fasta))
        #kcounter.count_kmers_custom(fasta, 20)
        #pass

    LOG.debug("Memory after kmer counting = {}G".format(get_memory_usage()))

    # spawn processes for amplicon finding before loading all the conserved kmers into memory 
    input_queue = multiprocessing.Queue()
    output_queue = multiprocessing.Queue(10000)
    threads = 7
    max_dist = 400
    min_dist = 100
    max_mem = 60
    pool = create_worker_pool(threads, input_queue, output_queue, k, max_dist, min_dist)
 

    # find conserved kmers
    conserved_kmers = []
    for kmer in kcounter.get_conserved_kmers(len(fastas)):
        conserved_kmers.append(kmer)

    LOG.info("Found {} fuzzy kmers conserved in all genomes.".format(len(conserved_kmers)))

    # send kcounter for garbage collection
    del kcounter


    #
    ## Generate all properly spaced (same contig, strand, and within some distance) conserved kmer pairs
    ## Put the properly spaced pairs in a data structure that doesn't care the order of the pairing (i.e. what strand the pair is found on per genome
    #

    amplicon_f = pooled_find_kmer_pairs(pool, conserved_kmers, input_queue, output_queue, threads, max_mem)

    #written_files = write_all_amplicons(Amplicon.existing_amplicons, fastas, k, "amplicon_fastas")
    written_files = write_amplicons_from_file(amplicon_f, fastas, k , "amplicon_fastas")

    get_amplicon_stats(written_files, fastas, k)


    """
    All that is left to do is to make a table of which genomes have unique sequences for each amplicon.

    I'm thinking there will be multiple files ideally. 

    1. Amplicon stats - sorted file that has # of Ns in fwd and rev primers and the number of genomes it can tell appart

    2. Amp X genome matrix - matrix that tells which genomes each amplicon can tell apart

    3. Fasta file for each amplicon that can tell apart any number of genomes uniquely (n > 1).
    """
    

def main_database(fastas, k, threads=8):
    counter = KmerCounterSQL("kmer_db", k, threads)

    counter.count_kmers(fastas)
    counter.generate_fuzzy_kmers()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-fastas", help="a bunch of fastas to compare", nargs="+", required=True)
    parser.add_argument("-k", help="value of k to use", required=True, type=int)
    args = parser.parse_args()

    #main(args.fastas, args.k)
    main_database(args.fastas, args.k)
