
import argparse
import os
import sys
import math
import time
import numpy as np
import logging
import re
import pickle
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
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


I also suspect that this program may not be accounting for reverse complements.

TODO: Print some final stats about amplicons (number unique, etc).

TODO: Look into speeding up queries
    
    - index all the join on columns -- this should work fairly well

    - reorder clauses to eliminate more potential results quicker

TODO: To speed up building the main database, I could not build the kmer table until the end and then select distinct. Foreign keys are not enabled by default so it wouldn't matter. 

NOTE: Changing the sqlite3 row factory to sqlite3.Row messes up multiprocessing. If I need the Row class, I should make sure to set it back to the default None before moving on. 


Performance Testing:
    Adding genomes to locations may have resulted in a very small (or no) performance increase but the size added to the database is negligible so it is worth it for just simplifying queries

"""

class Database(object):
    """ Serves as the base class for database objects """

    @staticmethod
    def open_database(database_f):
        """ Opens a database and returns a connection """
        return sqlite3.connect(database_f, check_same_thread=False)

    @staticmethod
    def init_database(conn):
        """ Set up some pragmas that will hopefully increase speed """
        conn.execute("PRAGMA main.page_size = 4096")
        conn.execute("PRAGMA main.cache_size = 5000")
        conn.execute("PRAGMA main.locking_mode = EXCLUSIVE")
        conn.execute("PRAGMA main.synchronous = NORMAL")
        #conn.execute("PRAGMA main.journal_mode = WAL")
        #conn.exeucte("PRAGMA main.temp_store = MEMORY")

    @staticmethod
    def _sql_results_to_file(fh, cursor):
        """ Writes results to a tab delimited file """
        while True:
            results = cursor.fetchmany(100000)
            if not results:
                break
            else:
                for result in results:
                    fh.write("\t".join([str(i) for i in result]) + "\n")

    @staticmethod
    def _iter_from_file(f, batchsize):
        """ Yields groups of at most batchsize from a file. Each group is a list of tuples """
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

    
    def lookup_param(self, key):
        """ Lookup a param in the Params table by a given key. Returns a list for keys that have multiple values, a single value for keys that have one value, or None for keys that are not in the database """

        results = self.dbc.execute("SELECT value FROM Params WHERE key = '{}'".format(key)).fetchall()

        if len(results) == 0:
            return None
        elif len(results) == 1:
            # unpack the first result (a tuple) and then unpack the first (and only) element from the tuple
            # try to return it as a type that makes sense (float for number, else string)
            try:
                return float(results[0][0])
            except ValueError:
                return results[0][0]
        elif len(results) > 1:
            # unpack each result and try to convert to sensible types
            all_results = []
            for result in results:
                try:
                    all_results.append(float(result[0]))
                except ValueError:
                    all_results.append(result[0])

            return all_results


class MainDatabase(Database):
    """
    Represents the main kmer database.

    Tables are Genomes, Contigs, Sequences, Locations, FuzzyKmers, and KmerToFuzzy
    """

    def __init__(self, database, k, threads, mismatches=1, stable_bp=2):
        self.database_f = database
        self.k = k
        self.threads = threads
        self.mismatches = mismatches

        # stables are # bases at beginning and end that cannot be mismatches
        self.stable_bp = stable_bp

        # open/make the database
        if os.path.isfile(self.database_f): 
            self.dbc = self.open_database(self.database_f)
        else:
            self.dbc = self.open_database(self.database_f)
            self.init_database(self.dbc)
            self.make_SQL_tables()


    def check_recount(self):
        """ Returns True if kmers need to be recounted """

        if self.k == self.lookup_param("k"):
            return False
        else:
            return True

    def check_refuzzify(self):
        """ Returns True is fuzzy kmers need to be recalculated """

        if self.k == self.lookup_param("k") and self.mismatches == self.lookup_param("mismatches") and self.stable_bp == self.lookup_param("stable_bp"):
            return False
        else:
            return True

    def make_SQL_tables(self):
        """ Makes SQL tables corresponding to all the info I want """

        # create table that holds the params
        self.dbc.execute(""" CREATE TABLE Params (
                        key TEXT,
                        value TEXT,
                        UNIQUE(key, value) ON CONFLICT IGNORE
                        );
                        """)

        self.dbc.executemany(""" INSERT INTO Params (key, value) VALUES (?, ?)""", (("k", self.k), ("mismatches", self.mismatches), ("stable_bp", self.stable_bp)))

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

        self.dbc.execute(""" CREATE TABLE Sequences (
                        id INTEGER PRIMARY KEY NOT NULL,
                        contig INTEGER NOT NULL,
                        seq TEXT,
                        FOREIGN KEY(contig) REFERENCES Contigs(id)
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
                        genome INTEGER NOT NULL,
                        contig INTEGER NOT NULL,
                        start INTEGER,
                        strand INTEGER,
                        kmer INTEGER NOT NULL,
                        FOREIGN KEY (genome) REFERENCES Genomes(id)
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

    def count_kmers(self, fastas, recalculate=False):
        """
        The counter here should do nothing but update the database.

        In multiple other threads the counter should be running counting kmers and should return multiple lists of tuples. The first contains all the unique kmers that need to be in kmers and the second contains Locations for the counted segment.
        """

        if self.check_recount():
            if recalculate:
                LOG.warning("Kmers need to be recounted. Got the -force option. Overwriting existing db.")
                self.dbc.execute("DELETE FROM Genomes")
                self.dbc.execute("DELETE FROM Contigs")
                self.dbc.execute("DELETE FROM Sequences")
                self.dbc.execute("DELETE FROM Kmers")

            else:
                raise ValueError("Kmers need to be recounted (current params don't match database params) but -force was not supplied.")

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


                    # add contig and sequence to the database and get the contig index
                    contig_id = self._add_contig(record.description, genome_id, str(record.seq))

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
            self.dbc.commit()

        for _ in workers:
            input_queue.put("STOP")

        for worker in workers:
            LOG.debug("Waiting for {} to join...".format(worker.name))
            worker.join()

    #
    ## Populating fuzzy kmer table
    #
    def generate_fuzzy_kmers(self, recalculate=False):
        """ Populates the FuzzyKmer and KmerToFuzzy tables 
        
        Checks first if any new genomes have been added and skips this step if none have.

        Fuzzy kmers can be recalculated using recalculate=True
        """

        refuzzify = self.check_refuzzify()

        # first check if fuzzy kmers need to be recalculated
        if refuzzify:
            if recalculate:
                LOG.info("Recalculating fuzzy kmers from scratch...")
                self.dbc.execute("DELETE FROM FuzzyKmers")
                self.dbc.execute("DELETE FROM KmerToFuzzy")
 
            else:
                raise ValueError("Fuzzy Kmers need to be recalculated(current params don't match database params) which would erase original data but -force was not supplied.")


        # this step takes a long time on huge tables so I want to skip it if possible
        # I think I want to bypass this in the interest of ensuring the fuzzy kmers get calculated when they should
        """
        if refuzzify is False and self.changed is False:
            LOG.info("No database changes detected. Skipping updating fuzzy kmers.")
            return
        """
            
        # select kmers that do not have an entry in the KmersToFuzzy table (ones that have not yet been processed)
        LOG.info("Selecting kmers to fuzzify...")
        cursor = self.dbc.execute("SELECT Kmers.id, Kmers.name from Kmers LEFT JOIN KmerToFuzzy ON kmers.id = KmerToFuzzy.kmer WHERE KmerToFuzzy.kmer is NULL")

        LOG.info("Fuzzifying and inserting fuzzy kmers into database...")
        result_num = 0
        for fuzzy_kmers, links in lazy_imap(self.threads, self._generate_fuzzy_entries, self._iter_sql_results(cursor, 500000)):

            result_num += 1
            LOG.debug("After getting result {}".format(result_num))
            get_memory_usage(verbose=True)

            LOG.debug("Inserting {} elements into FuzzyKmers...".format(len(fuzzy_kmers)))
            self.dbc.executemany("INSERT OR IGNORE INTO FuzzyKmers (id, name) VALUES (?, ?)", fuzzy_kmers)
            LOG.debug("Inserting {} elements into KmerToFuzzy...".format(len(links)))
            self.dbc.executemany("INSERT INTO KmerToFuzzy (kmer, fuzzy) VALUES (?, ?)", links)

        self.dbc.commit()
        
    def _generate_fuzzy_entries(self, kmers):
        """ Takes an iterable of kmers (id, name) and returns a list of fuzzy kmers ready to be inserted into the FuzzyKmer table and a list of links between FuzzyKmers and real kmers"""

        #LOG.debug("Beginning fuzzifying batch...")
        fuzzy_kmers = []
        links = []
        for kmer_id, kmer in kmers:
            for fuzzy_kmer in self._fuzzify_kmer(kmer):
                fuzzy_id_str = fuzzy_kmer.replace("A", "0").replace("C", "1").replace("G", "2").replace("T", "3").replace("N", "4")
                fuzzy_id = int(fuzzy_id_str, 5)
            
                fuzzy_kmers.append((fuzzy_id, fuzzy_kmer))
                links.append((kmer_id, fuzzy_id))

        #LOG.debug("Done fuzzifying batch...")
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

        cursor = self.dbc.execute("INSERT INTO Genomes (name) VALUES (?)", (genome_name,))
        return cursor.lastrowid

    def _add_contig(self, contig_name, genome_id, sequence):
        """ Adds a contig and its sequence and returns the contig id in the database """
        
        # add the contig
        cursor = self.dbc.execute("INSERT INTO Contigs (name, genome) VALUES (?, ?)", (contig_name, genome_id))
        contig_id = cursor.lastrowid

        # add the sequence
        self.dbc.execute("INSERT INTO Sequences (contig, seq) VALUES (?, ?)", (contig_id, sequence))

        return contig_id

    def _add_kmers(self, kmers):
        """ Insert ignore new kmers from a list of kmer iterables """
        self.dbc.executemany('INSERT OR IGNORE INTO Kmers (id, name) VALUES (?, ?)', kmers)

    def _add_locations(self, locations):
        """ Insert new locations from a list of location tuples """
        self.dbc.executemany('INSERT INTO Locations (genome, contig, start, strand, kmer) VALUES (?, ?, ?, ?, ?)', locations)


    #
    ## Run-specific commands 
    #

    def prepare_run_db(self):
        self.open_database

    def find_conserved_kmers(self, genomes=None, fuzzy=False):
        """ Finds conserved kmers in all (default) or a specific set of genomes (specified by genomes argument) 
        It would be nice to somehow filter out duplicates where duplicates = :
            1. kmers that are consecutive in all genomes
            2. fuzzy kmers that are logical duplicates (each perfectly matching kmer has something like 16 logical duplicates).
        """
        
        # make a temporary table with only the genomes we want to include
        # If genomes not specified, this table will be an exact clone of Genomes
        self.dbc.execute(""" DROP TABLE IF EXISTS SubGen """)
        self.dbc.execute("""CREATE TEMPORARY TABLE SubGen (
                        id INTEGER PRIMARY KEY NOT NULL,
                        name STRING
                        );
                        """)
        if genomes:
            self.dbc.execute("""INSERT INTO SubGen (id, name)
                            SELECT id, name FROM Genomes
                            WHERE Genomes.id IN ({ph})
                            """.format(ph=",".join(["?"]*len(genomes))), genomes)
        else:
            self.dbc.execute("""INSERT INTO SubGen (id, name)
                            SELECT id, name FROM Genomes
                            """)
 
        return

        # Make a temporary table to hold the conserved Kmers
        self.dbc.execute("DROP TABLE IF EXISTS ConservedKmers")

        self.dbc.execute("""
                CREATE TABLE ConservedKmers (
                id INTEGER PRIMARY KEY NOT NULL,
                kmer INTEGER NOT NULL,
                genome INTEGER NOT NULL,
                contig INTEGER NOT NULL,
                start INTEGER NOT NULL,
                strand INTEGER NOT NULL
                );
                """)


        if fuzzy:
            # select fuzzy kmers that are present in all genomes
            LOG.debug("Finding conserved fuzzy kmers...")

            self.dbc.execute("""
                    -- Insert results into the temp table
                    INSERT INTO ConservedKmers (kmer, genome, contig, start, strand) 
                    
                    -- Select from Location
                    SELECT cons_kmer, Contigs.genome, Locations.contig, Locations.start, Locations.strand 
                    FROM Locations 

                    -- join Contigs to be able to join SubGen
                    INNER JOIN Contigs 
                        ON Contigs.id = Locations.contig 
                    
                    -- join SubGen to check genome number
                    INNER JOIN SubGen
                        ON Contigs.genome = SubGen.id


                    -- join with a select that checks  all genomes represented
                    INNER JOIN (
                            SELECT FuzzyKmers.id cons_kmer, KmerToFuzzy.kmer as fkmer
                            FROM FuzzyKmers

                            -- join KmerToFuzzy to get the kmers each fuzzy represents
                            INNER JOIN KmerToFuzzy
                                ON FuzzyKmers.id = KmerToFuzzy.fuzzy

                            -- join locations to get counts
                            INNER JOIN Locations
                                ON KmerToFuzzy.kmer = Locations.kmer

                            -- join contigs to get genomes
                            INNER JOIN Contigs 
                                ON Locations.contig = Contigs.id 

                            -- group by fuzzy kmer
                            GROUP BY FuzzyKmers.id

                            -- contrain to all genomes represented
                            HAVING COUNT(DISTINCT Contigs.genome) = (SELECT COUNT(*) FROM SubGen)) 
                        ON Locations.kmer = fkmer

               """)


        else:

            # select kmers that are present in all genomes
            LOG.debug("Finding conserved kmers...")

            self.dbc.execute("""
                    -- insert results into temp table
                    INSERT INTO ConservedKmers (kmer, genome, contig, start, strand) 
                    
                    -- select from locations
                    SELECT Locations.kmer, Contigs.genome, Locations.contig, Locations.start, Locations.strand 
                    FROM Locations 

                    -- join contigs to be able to join SubGen
                    INNER JOIN Contigs 
                        ON Contigs.id = Locations.contig 

                    -- join subgen to limit to only required genomes
                    INNER JOIN SubGen
                        ON Contigs.genome = SubGen.id

                    -- join to select that check for kmers in all genomes
                    INNER JOIN (
                            SELECT Locations.kmer ckmer 
                            FROM Locations 

                            -- join contigs to get genome
                            JOIN Contigs 
                                ON Locations.contig = Contigs.id

                            -- collapse around kmers
                            GROUP BY Locations.kmer 

                            -- contrain to only kmers found in all genomes
                            HAVING COUNT(DISTINCT Contigs.genome) = (SELECT COUNT(*) FROM SubGen)) 
                        ON Locations.kmer = ckmer
                """)

        self.database.commit()
        #cursor = self.dbc.execute("SELECT * FROM ConservedKmers")
        #print(cursor.fetchall())

    def find_potential_amplicons(self, kmer_table="Kmers", min_genomes=1, min_dist=100, max_dist=400):
        """ Finds all potential amplicons from conserved kmers """

        LOG.info("Making amplicons table...")

        # -- create temporary table 
        self.dbc.execute("""
                -- make temporary table
                CREATE TEMPORARY TABLE temp_amplicons AS 
                
                -- populate it from this select
                SELECT A.kmer AS kmer1, B.kmer AS kmer2, A.start as start1, B.start as start2, A.genome AS genome, A.contig AS contig, A.strand AS strand 
                
                FROM ConservedKmers A 
                
                -- self join to compare pairwise
                INNER JOIN ConservedKmers B 
                    ON (
                        -- contig and strand must be the same
                        A.contig = B.contig AND 
                        A.strand = B.strand AND 
                        
                        -- distance must be appropriate and consistent with strand
                        (
                            (A.strand = '1' AND B.start - A.start BETWEEN {min_dist} AND {max_dist}) 
                        OR 
                            (A.strand = '-1' AND A.start - B.start BETWEEN {min_dist} AND {max_dist})
                        )

                        -- no self matches
                        AND A.kmer != B.kmer
                    )
                """.format(min_dist=min_dist, max_dist=max_dist))

        LOG.debug("temp_amplicons made... finishing table...")


        # -- query temporary table
        self.dbc.execute("DROP TABLE IF EXISTS Amplicons")
        # this statement appends columns for num_genomes and num_total
        self.dbc.execute("""
                -- create a table
                CREATE TABLE Amplicons AS 
                
                    -- populate it with the results from this select
                    SELECT K1.name AS kmer1, K2.name AS kmer2, tmp_amp.start1, tmp_amp.start2, tmp_amp.genome, tmp_amp.contig, tmp_amp.strand, num_genomes, num_total 
                    
                    FROM temp_amplicons tmp_amp 
                   
                    -- join with some counts
                    INNER JOIN 
                        (SELECT temp_amplicons.kmer1 AS ta1, temp_amplicons.kmer2 AS ta2, COUNT(DISTINCT temp_amplicons.genome) AS num_genomes, COUNT(*) AS num_total 
                    
                        FROM temp_amplicons 
                        
                        -- group by amplicon
                        GROUP BY temp_amplicons.kmer1, temp_amplicons.kmer2 
                        
                        -- require a minimum number of genomes
                        HAVING num_genomes >= {min_genomes}) 
                    
                        -- merge counts with original table
                        ON tmp_amp.kmer1 = ta1 AND tmp_amp.kmer2 = ta2 
                    
                    -- join with specified kmer table to get kmer1 sequence
                    INNER JOIN {kmer_table} as K1 
                        ON tmp_amp.kmer1 = K1.id 
                    
                    -- join with specified kmer table to get kmer2 sequence
                    INNER JOIN {kmer_table} as K2 
                        ON tmp_amp.kmer2 = K2.id 
                        """.format(kmer_table=kmer_table, min_genomes=min_genomes))

        # drop the temporary table
        self.dbc.execute("DROP TABLE temp_amplicons")

        # add an index to amplicons
        self.dbc.execute("CREATE INDEX amplicon_indx on Amplicons (kmer1, kmer2)")

        self.database.commit()

    def get_amplicon_stats(self):

        LOG.info("Checking amplicons and getting stats...")
        cursor = self.dbc.execute("SELECT SubGen.name FROM SubGen ORDER BY SubGen.name")
        genome_names = [str(itm[0]) for itm in cursor.fetchall()]

      
        # set some stats for writing amplicons
        self.amp_write_genomes = len(genome_names)
        self.amp_write_unique = len(genome_names)
        self.amp_write_duplicates = 0
        self.amp_write_N = 99


        # make amplicon directory if it doesn't exist
        if not os.path.isdir(self.amplicon_dir):
            os.mkdir(self.amplicon_dir)

        ordered_matches = "SELECT Amplicons.kmer1, Amplicons.kmer2 FROM Amplicons GROUP BY Amplicons.kmer1, Amplicons.kmer2 ORDER BY num_genomes DESC, num_total ASC"

        cursor = self.dbc.execute(ordered_matches)

        amplicon_stats = {}
        for result in lazy_imap(self.threads, self.check_amplicon, self.generate_amp_packages(cursor)):
            print(result)
            amplicon_stats.update(result)


        # sort the amplicons by num unique_genomes number (descending) and number of N's(ascending), and then num duplicates
        sorted_amplicons = sorted(amplicon_stats, key=lambda amp_id: (-1 * amplicon_stats[amp_id]["num_unique"], amplicon_stats[amp_id]["fwd_N"] + amplicon_stats[amp_id]["rev_N"], amplicon_stats[amp_id]["num_duplicates"]))


        with open(self.amplicon_stats_f, 'w') as STATS, open(self.amplicon_matr_f, 'w') as MATR:
            # write headers to both files
            STATS.write("\t".join(("amplicon", "num_genomes", "num_unique", "num_duplicates", "fwd_N", "rev_N")) + "\n")
            MATR.write("\t".join(["amplicon"] + list(genome_names)) + "\n")

            for amp in sorted_amplicons:
                # k for key not length of kmer
                write_stats = [str(amp)] + [str(amplicon_stats[amp][k]) for k in ["num_genomes", "num_unique", "num_duplicates", "fwd_N", "rev_N"]]
                write_matr = [str(amp)]
                for name in genome_names:
                    try:
                        write_matr.append(";".join(amplicon_stats[amp]["genome_clusters"][name]))
                    except KeyError:
                        write_matr.append("0")

                STATS.write("\t".join(write_stats) + "\n")
                MATR.write("\t".join(write_matr) + "\n")

    def generate_amp_packages(self, cursor):
        """ Generates so-called "amplicon packages" from a cursor that selects amplicons """

        amp_id = 0
        for batch in self._iter_sql_results(cursor, 100):
            for amp in batch:
                amp_id += 1
                k1, k2 = amp

                seqs = self.dbc.execute("""
                        SELECT Genomes.name, Contigs.name, Amplicons.start1 as start, Amplicons.start2 + {k} as end, Amplicons.strand, SUBSTR(Sequences.seq, Amplicons.start1 + 1, Amplicons.start2 + {k} - Amplicons.start1) 
                        
                        FROM Amplicons 
                        
                        INNER JOIN Genomes 
                            ON Amplicons.genome = Genomes.id 
                        
                        INNER JOIN Contigs 
                            ON Amplicons.contig = Contigs.id 
                        
                        INNER JOIN Sequences 
                            ON Amplicons.contig = Sequences.id WHERE Amplicons.kmer1 = '{k1}' 
                            AND Amplicons.kmer2 = '{k2}'
                        
                        """.format(k=self.k, k1=k1, k2=k2)).fetchall()
                yield (amp_id, k1, k2, seqs)

    def check_amplicon(self, amp_package):
        """ Checks amplicons and returns a dict of amplicon stats """
        
        # unpack the package
        amp_id, fwd_prim, rev_prim, db_sequences = amp_package

        amplicon_name = "amp{amp_id}_{fwd_prim}-{rev_prim}.fasta".format(amp_id=amp_id, fwd_prim=fwd_prim, rev_prim=rev_prim)
        # initialize the results dict
        amplicon_stats = {amplicon_name: {}}

        amplicon_stats[amplicon_name]["fwd_N"] = fwd_prim.count("N")
        amplicon_stats[amplicon_name]["rev_N"] = rev_prim.count("N")

        #
        ## convert database data to BioPython Sequence objects
        #
        full_seqs = []
        inter_primer_only = {}
        for db_data in db_sequences:
            genome, contig, start, end, strand, seq = db_data
            
            header = "genome={genome};contig={contig};location=[{start}:{end}];strand={strand}".format(genome=genome, contig=contig, start=start, end=end, strand=strand)
            seq_obj = SeqRecord(seq=Seq(seq), id=header, description="")

            # reverse complement if necessary
            if strand == "-1":
                seq_obj = seq_obj.reverse_complement()


            # store one copy in full_seqs and one copy in inter_primer_only
            full_seqs.append(seq_obj)
            inter_primer_only[(genome, header)] = str(seq_obj[self.k:-self.k].seq)
 
        #
        ## Place each sequence into a cluster based on uniqueness
        #
        clusters = {}
        cluster_counts = {}
        cluster = 0
        for key1, seq1 in inter_primer_only.items():

            # skip seqs that already have cluster assigned
            if key1 in clusters:
                continue
            else:
                cluster += 1
                clusters[key1] = str(cluster)
                cluster_counts[cluster] = 1

            for key2, seq2 in inter_primer_only.items():
                if key1 is key2:
                    continue

                # skip seqs that have already been processed
                elif key2 in clusters:
                    continue

                else:
                    # check if the seqs are == 
                    if seq1 == seq2:
                        clusters[key2] = str(cluster)
                        cluster_counts[cluster] += 1

        # Find genomes that are parts of unique clusters
        duplicated_genomes = 0
        unique_genomes = set()     
        genome_clusters = {}
        for key in clusters:
            genome, header = key

            cluster = clusters[key]

            if cluster_counts[int(cluster)] == 1:
                unique_genomes.add(genome)
            try:
                genome_clusters[genome].append(clusters[key])
                duplicated_genomes += 1
            except KeyError:
                genome_clusters[genome] = [clusters[key]]

        amplicon_stats[amplicon_name]["num_genomes"] = len(genome_clusters)
        amplicon_stats[amplicon_name]["num_unique"] = len(unique_genomes)
        amplicon_stats[amplicon_name]["num_duplicates"] = duplicated_genomes
        amplicon_stats[amplicon_name]["genome_clusters"] = genome_clusters

        # check some standards to determine if this amp should be written
        if amplicon_stats[amplicon_name]["num_genomes"] >= self.amp_write_genomes and \
                amplicon_stats[amplicon_name]["num_unique"] >= self.amp_write_unique and \
                amplicon_stats[amplicon_name]["num_duplicates"] <= self.amp_write_duplicates and \
                amplicon_stats[amplicon_name]["fwd_N"] + amplicon_stats[amplicon_name]["rev_N"] <= self.amp_write_N:

            with open(self.amplicon_dir + "/" + amplicon_name, 'w') as OUT:
                SeqIO.write(full_seqs, OUT, "fasta")

        return amplicon_stats


class DatabaseRun(Database):
    """ Represents a single query/run of the database 
    
    I need a new way of getting/setting/saving params. I want the run_db.Params table to represent everything about the run:

        - genomes included
        - fuzzy
        - amp min len
        - amp max len
        - amp frac genomes

        and perhaps also the write params.

    Given that genomes included will be more than one, I think I need to change the Params schema to be a k: v system where I have something like index:key:value and a table would look like this:

        1:genome:Gen1
        2:genome:Gen2
        3:fuzzy:True
        4:amp_min_len:100
        5:amp_max_len:400
        etc...
    
    """
    
    def __init__(self, main_db_f, output_dir, run_prefix, threads): 

        self.threads = threads

        #
        ## path vars
        #
        self.output_dir = output_dir
        self.prefix = run_prefix
        self.run_database_f = os.path.join(self.output_dir, self.prefix + ".db")
        self.amplicon_dir = os.path.join(self.output_dir, self.prefix + ".amplicons", "")
        self.amplicon_stats_f = os.path.join(self.output_dir, self.prefix + ".amplicon_stats.txt")
        self.amplicon_matr_f = os.path.join(self.output_dir, self.prefix + ".amplicon_matrix.txt")


        self.dbc = self.open_database(main_db_f)

    def _attach_run_db(self, force):
        """ Attached the run_specific db and checks if the file exists."""

        if os.path.exists(self.run_database_f):
            # if force, remove the run_specific database
            if force:
                LOG.info("Force option received: Overwriting run_specific tables if they exist.")
                os.unlink(self.run_database_f)
            else:
                raise ValueError("Run databse already exists. Refusing to overwrite without the -force option. It is recommended to just specify a different -prefix")

        # attach the run_specific database
        self.dbc.execute("ATTACH '{}' AS run_specific".format(self.run_database_f))
 
        # make params table if it doesn't exist
        exists = self.dbc.execute("SELECT name FROM run_specific.sqlite_master WHERE type='table' AND name='Params'").fetchone()

        if exists is None:
            # run is new, make the params table
            self.dbc.execute("CREATE TABLE run_specific.Params (key, value, UNIQUE (key, value) ON CONFLICT IGNORE)")
            self.dbc.commit()

    def add_params_to_database(self, **params):
        """ adds params to database where params is a bunch of keyword args """
        self.dbc.executemany("INSERT INTO run_specific.params (key, value) VALUES (?, ?)", [(k, v) for k, v in params.items()])
        
        self.dbc.commit()

    def lookup_param(self, key):
        """ Lookup a param in the Params table by a given key. Returns a list for keys that have multiple values, a single value for keys that have one value, or None for keys that are not in the database """

        results = self.dbc.execute("SELECT value FROM run_specific.Params WHERE key = '{}'".format(key)).fetchall()

        if len(results) == 0:
            return None
        elif len(results) == 1:
            # unpack the first result (a tuple) and then unpack the first (and only) element from the tuple
            # try to return it as a type that makes sense (float for number, else string)
            try:
                return float(results[0][0])
            except ValueError:
                return results[0][0]
        elif len(results) > 1:
            # unpack each result and try to convert to sensible types
            all_results = []
            for result in results:
                try:
                    all_results.append(float(result[0]))
                except ValueError:
                    all_results.append(result[0])

            return all_results

    def check_previous_run(self):
        """ Checks the params used in a previous run in the same run_db to see if it is possible to reuse some data. This is not implemented yet.
        """
        # before starting, it is important to check a few basic things to see if any of the run can be resused

        # compare run k to main_db k
        if super().lookup_param("k") != self.lookup_param("k"):
            LOG.warning("Main database was created using a different k than this run database.")
            
        # compare run genomes to current genomes
        if set(self.lookup_param("genomes")) != set(genomes):
            LOG.warning("Genomes specified for this run are different than genomes specified for the previous run done in this database.")
        
        # compare run fuzzy to current fuzzy

        # compare run mismatches to main mismatches

        # compare run stable_bp to main stable_bp

    def process(self, genomes, fuzzy, amp_min_len, amp_max_len, amp_frac_genomes, force):
        """ Process the run recalculating as few steps as possible """

        #self._attach_run_db(force)
        self.dbc.execute("ATTACH '{}' as run_specific".format(self.run_database_f))

        # set genomes = to all genomes in the database if none were specified
        if genomes is None:
            cursor = self.dbc.execute("SELECT name FROM Genomes")
            genomes = []
            for genome in cursor.fetchall():
                genomes.append(genome[0])


        # add the params to the database
        self.add_params_to_database(fuzzy=fuzzy, amp_min_len=amp_min_len, amp_max_len=amp_max_len, amp_frac_genomes=amp_frac_genomes)

        for genome in genomes:
            self.add_params_to_database(genomes=genome)


        # analyze the table to help create the best plan
        #self.dbc.execute("ANALYZE")

        #self.find_conserved_kmers()
        self.find_potential_amplicons()

        self.get_amplicon_stats()

    def find_conserved_kmers(self):
        """ Finds conserved kmers in all (default) or a specific set of genomes (specified by genomes argument) 
        It would be nice to somehow filter out duplicates where duplicates = :
            1. kmers that are consecutive in all genomes
            2. fuzzy kmers that are logical duplicates (each perfectly matching kmer has something like 16 logical duplicates).
        """
        
        # make a table with only the genomes we want to include
        self.dbc.execute("""CREATE TABLE run_specific.SubGen (
                        id INTEGER PRIMARY KEY NOT NULL,
                        name STRING
                        );
                        """)

        # lookup the genomes to include from the params table
        genomes = self.lookup_param("genomes")
        # convert to a list if not already a list
        if type(genomes) != list:
            genomes = [genomes]

        command = """
                        INSERT INTO SubGen (id, name)
                        SELECT id, name FROM Genomes
                        WHERE Genomes.name IN ('{placeholder}')
                        """.format(placeholder="', '".join([genome for genome in genomes]))
 

        cursor = self.dbc.execute("""
                        INSERT INTO SubGen (id, name)
                        SELECT id, name FROM Genomes
                        WHERE Genomes.name IN ({placeholder})
                        """.format(placeholder=",".join(["?"]*len(genomes))), genomes)
 
        # Make a table to hold the conserved Kmers
        self.dbc.execute("""
                CREATE TABLE run_specific.ConservedKmers (
                id INTEGER PRIMARY KEY NOT NULL,
                kmer INTEGER NOT NULL,
                genome INTEGER NOT NULL,
                contig INTEGER NOT NULL,
                start INTEGER NOT NULL,
                strand INTEGER NOT NULL
                );
                """)

        # lookup the fuzzy value from the params table
        fuzzy = self.lookup_param("fuzzy")


        if fuzzy:
            # select fuzzy kmers that are present in all genomes
            LOG.debug("Finding conserved fuzzy kmers...")

            cursor = self.dbc.execute("""
                    -- Insert results into the temp table
                    INSERT INTO ConservedKmers (kmer, genome, contig, start, strand) 

                    -- Select from Location
                    SELECT cons_kmer, Locations.genome, Locations.contig, Locations.start, Locations.strand 
                    FROM Locations 

                    -- join SubGen to check genome number
                    INNER JOIN SubGen
                        ON Locations.genome = SubGen.id


                    -- join with a select that checks  all genomes represented
                    INNER JOIN (
                            SELECT FuzzyKmers.id cons_kmer, KmerToFuzzy.kmer as fkmer
                            FROM FuzzyKmers

                            -- join KmerToFuzzy to get the kmers each fuzzy represents
                            INNER JOIN KmerToFuzzy
                                ON FuzzyKmers.id = KmerToFuzzy.fuzzy

                            -- join locations to get counts
                            INNER JOIN Locations
                                ON KmerToFuzzy.kmer = Locations.kmer

                            -- group by fuzzy kmer
                            GROUP BY FuzzyKmers.id

                            -- contrain to all genomes represented
                            HAVING COUNT(DISTINCT Locations.genome) = (SELECT COUNT(*) FROM SubGen)) 
                        ON Locations.kmer = fkmer

               """)
        else:

            # select kmers that are present in all genomes
            LOG.debug("Finding conserved kmers...")

            self.dbc.execute("""
                    -- insert results into temp table
                    INSERT INTO ConservedKmers (kmer, genome, contig, start, strand) 
                    
                    -- select from locations
                    SELECT Locations.kmer, Contigs.genome, Locations.contig, Locations.start, Locations.strand 
                    FROM Locations 

                    -- join contigs to be able to join SubGen
                    INNER JOIN Contigs 
                        ON Contigs.id = Locations.contig 

                    -- join subgen to limit to only required genomes
                    INNER JOIN SubGen
                        ON Contigs.genome = SubGen.id

                    -- join to select that check for kmers in all genomes
                    INNER JOIN (
                            SELECT Locations.kmer ckmer 
                            FROM Locations 

                            -- join contigs to get genome
                            JOIN Contigs 
                                ON Locations.contig = Contigs.id

                            -- collapse around kmers
                            GROUP BY Locations.kmer 

                            -- contrain to only kmers found in all genomes
                            HAVING COUNT(DISTINCT Contigs.genome) = (SELECT COUNT(*) FROM SubGen)) 
                        ON Locations.kmer = ckmer
                """)

        self.dbc.commit()

    def find_potential_amplicons(self):
        """ Finds all potential amplicons from conserved kmers 
        
        TODO: Add k to length calculations to get actual length
        """

        # lookup fuzzy so we know which table to use
        fuzzy = self.lookup_param("fuzzy")
        if fuzzy:
            kmer_table = "FuzzyKmers"
        else:
            kmer_table = "Kmers"

        LOG.info("Making amplicons table...")


        amp_min_len = self.lookup_param("amp_min_len")
        amp_max_len = self.lookup_param("amp_max_len")


        # I might could speed this up by breaking up the OR into two separate inserts

        # -- create temporary table 
        self.dbc.execute("""
                -- make temporary table
                CREATE TABLE run_specific.temp_amplicons AS 
                
                -- populate it from this select
                SELECT A.kmer AS kmer1, B.kmer AS kmer2, A.start as start1, B.start as start2, A.genome AS genome, A.contig AS contig, A.strand AS strand 
                
                FROM ConservedKmers A 
                
                -- self join to compare pairwise
                INNER JOIN ConservedKmers B 
                    ON (
                        -- contig and strand must be the same
                        -- contig first because it is more specific than strand
                        A.contig = B.contig AND 
                        A.strand = B.strand AND 

                        -- no self matches
                        A.kmer != B.kmer AND
                        
                        -- distance must be appropriate and consistent with strand
                        (
                            (A.strand = '1' AND B.start - A.start BETWEEN {amp_min_len} AND {amp_max_len}) 
                        OR 
                            (A.strand = '-1' AND A.start - B.start BETWEEN {amp_min_len} AND {amp_max_len})
                        )

                   )
                """.format(amp_min_len=amp_min_len, amp_max_len=amp_max_len))

        LOG.debug("temp_amplicons made... finishing table...")


        # get the number of genomes that must be present in the amplicon
        total_genomes = self.dbc.execute("SELECT COUNT(*) FROM run_specific.SubGen").fetchone()[0]
        amp_frac_genomes = self.lookup_param("amp_frac_genomes")
        min_num_genomes = math.ceil(total_genomes / amp_frac_genomes)

        # this statement appends columns for num_genomes and num_total
        self.dbc.execute("""
                -- create a table
                CREATE TABLE run_specific.Amplicons AS 
                
                    -- populate it with the results from this select
                    SELECT K1.name AS kmer1, K2.name AS kmer2, tmp_amp.start1, tmp_amp.start2, tmp_amp.genome, tmp_amp.contig, tmp_amp.strand, num_genomes, num_total 
                    
                    FROM temp_amplicons tmp_amp 
                   
                    -- join with some counts
                    INNER JOIN 
                        (SELECT temp_amplicons.kmer1 AS ta1, temp_amplicons.kmer2 AS ta2, COUNT(DISTINCT temp_amplicons.genome) AS num_genomes, COUNT(*) AS num_total 
                    
                        FROM temp_amplicons 
                        
                        -- group by amplicon
                        GROUP BY temp_amplicons.kmer1, temp_amplicons.kmer2 
                        
                        -- require a minimum number of genomes
                        HAVING num_genomes >= {min_num_genomes}) 
                    
                        -- merge counts with original table
                        ON tmp_amp.kmer1 = ta1 AND tmp_amp.kmer2 = ta2 
                    
                    -- join with specified kmer table to get kmer1 sequence
                    INNER JOIN {kmer_table} as K1 
                        ON tmp_amp.kmer1 = K1.id 
                    
                    -- join with specified kmer table to get kmer2 sequence
                    INNER JOIN {kmer_table} as K2 
                        ON tmp_amp.kmer2 = K2.id 
                        """.format(kmer_table=kmer_table, min_num_genomes=min_num_genomes))

        # drop the temporary table
        #self.dbc.execute("DROP TABLE temp_amplicons")

        # add an index to amplicons
        self.dbc.execute("CREATE INDEX run_specific.amplicon_indx on Amplicons (kmer1, kmer2)")

        self.dbc.commit()

    def get_amplicon_stats(self):

        LOG.info("Checking amplicons and getting stats...")
        cursor = self.dbc.execute("SELECT SubGen.name FROM SubGen ORDER BY SubGen.name")
        genome_names = [str(itm[0]) for itm in cursor.fetchall()]

      
        ordered_matches = "SELECT rowid, Amplicons.kmer1, Amplicons.kmer2 FROM Amplicons GROUP BY Amplicons.kmer1, Amplicons.kmer2 ORDER BY num_genomes DESC, num_total ASC"

        cursor = self.dbc.execute(ordered_matches)

        amplicon_stats = {}
        for result in lazy_imap(self.threads, self.check_amplicon, self.generate_amp_packages(cursor)):
            amplicon_stats.update(result)


        # sort the amplicons by num unique_genomes number (descending) and number of N's(ascending), and then num duplicates
        sorted_amplicons = sorted(amplicon_stats, key=lambda amp_id: (-1 * amplicon_stats[amp_id]["num_unique"], amplicon_stats[amp_id]["fwd_N"] + amplicon_stats[amp_id]["rev_N"], amplicon_stats[amp_id]["num_duplicates"]))


        # make amplicon directory if it doesn't exist
        if not os.path.isdir(self.amplicon_dir):
            os.mkdir(self.amplicon_dir)


        with open(self.amplicon_stats_f, 'w') as STATS, open(self.amplicon_matr_f, 'w') as MATR:
            # write headers to both files
            STATS.write("\t".join(("amplicon", "num_genomes", "num_unique", "num_duplicates", "fwd_N", "rev_N")) + "\n")
            MATR.write("\t".join(["amplicon"] + list(genome_names)) + "\n")

            for indx, amp in enumerate(sorted_amplicons):
                
                # write the top 100 amplicons as fastas
                if indx < 100:
                    # trim the amp name to just the amp_id
                    amp_id = amp.split("_")[0][3:]

                    self.write_amplicon_fasta(amp_id, amp + ".fasta", self.amplicon_dir)


                # k for key not length of kmer
                write_stats = [str(amp)] + [str(amplicon_stats[amp][k]) for k in ["num_genomes", "num_unique", "num_duplicates", "fwd_N", "rev_N"]]
                write_matr = [str(amp)]
                for name in genome_names:
                    try:
                        write_matr.append(";".join(amplicon_stats[amp]["genome_clusters"][name]))
                    except KeyError:
                        write_matr.append("0")

                STATS.write("\t".join(write_stats) + "\n")
                MATR.write("\t".join(write_matr) + "\n")
        
        LOG.info("Successfully finished generating amplicon stats.") 
        LOG.info("Wrote amplicon stats to...  {}".format(self.amplicon_stats_f))
        LOG.info("Wrote amplicon uniqueness matrix to...  {}".format(self.amplicon_matr_f))
        LOG.info("Wrote top amplicon FASTA files to...  {}".format(self.amplicon_dir))

    def generate_amp_packages(self, cursor):
        """ Generates so-called "amplicon packages" from a cursor that selects amplicons """

        for batch in self._iter_sql_results(cursor, 100):
            for amp in batch:
                amp_id, k1, k2 = amp
                seqs = self.get_amplicon_seqs(amp_id)

                yield (amp_id, k1, k2, seqs)

    def check_amplicon(self, amp_package):
        """ Checks amplicons and returns a dict of amplicon stats """
        
        k = int(super().lookup_param("k"))

        # unpack the package
        amp_id, fwd_prim, rev_prim, db_sequences = amp_package

        amplicon_name = "amp{amp_id}_{fwd_prim}-{rev_prim}".format(amp_id=amp_id, fwd_prim=fwd_prim, rev_prim=rev_prim)
        # initialize the results dict
        amplicon_stats = {amplicon_name: {}}

        amplicon_stats[amplicon_name]["fwd_N"] = fwd_prim.count("N")
        amplicon_stats[amplicon_name]["rev_N"] = rev_prim.count("N")

        #
        ## convert database data to BioPython Sequence objects
        #
        full_seqs = []
        inter_primer_only = {}
        for db_data in db_sequences:
            genome, contig, start, end, strand, seq = db_data
            
            header = "genome={genome};contig={contig};location=[{start}:{end}];strand={strand}".format(genome=genome, contig=contig, start=start, end=end, strand=strand)
            seq_obj = SeqRecord(seq=Seq(seq), id=header, description="")

            # reverse complement if necessary
            if strand == "-1":
                seq_obj = seq_obj.reverse_complement()


            # store one copy in full_seqs and one copy in inter_primer_only
            full_seqs.append(seq_obj)
            inter_primer_only[(genome, header)] = str(seq_obj[k:-k].seq)
 
        #
        ## Place each sequence into a cluster based on uniqueness
        #
        clusters = {}
        cluster_counts = {}
        cluster = 0
        for key1, seq1 in inter_primer_only.items():

            # skip seqs that already have cluster assigned
            if key1 in clusters:
                continue
            else:
                cluster += 1
                clusters[key1] = str(cluster)
                cluster_counts[cluster] = 1

            for key2, seq2 in inter_primer_only.items():
                if key1 is key2:
                    continue

                # skip seqs that have already been processed
                elif key2 in clusters:
                    continue

                else:
                    # check if the seqs are == 
                    if seq1 == seq2:
                        clusters[key2] = str(cluster)
                        cluster_counts[cluster] += 1

        # Find genomes that are parts of unique clusters
        duplicated_genomes = 0
        unique_genomes = set()     
        genome_clusters = {}
        for key in clusters:
            genome, header = key

            cluster = clusters[key]

            if cluster_counts[int(cluster)] == 1:
                unique_genomes.add(genome)
            try:
                genome_clusters[genome].append(clusters[key])
                duplicated_genomes += 1
            except KeyError:
                genome_clusters[genome] = [clusters[key]]

        amplicon_stats[amplicon_name]["num_genomes"] = len(genome_clusters)
        amplicon_stats[amplicon_name]["num_unique"] = len(unique_genomes)
        amplicon_stats[amplicon_name]["num_duplicates"] = duplicated_genomes
        amplicon_stats[amplicon_name]["genome_clusters"] = genome_clusters

        return amplicon_stats

    def write_amplicon_fasta(self, amplicon_id, file_name, output_dir):
        """ Writes an amplicon fasta in the specified directory """
        out_path = os.path.join(output_dir, file_name)

        seqs = self.get_amplicon_seqs(amplicon_id)
        seq_objs = []
        for seq_data in seqs:

            genome, contig, start, end, strand, seq = seq_data
            header = "genome={genome};contig={contig};location=[{start}:{end}];strand={strand}".format(genome=genome, contig=contig, start=start, end=end, strand=strand)
            seq_obj = SeqRecord(seq=Seq(seq), id=header, description="")

            # reverse complement if necessary
            if strand == "-1":
                seq_obj = seq_obj.reverse_complement()

            # add to the seq_objs list
            seq_objs.append(seq_obj)

        with open(out_path, 'w') as OUT:
            SeqIO.write(seq_objs, OUT, "fasta")

    def get_amplicon_seqs(self, amp_id):
        """ Returns a list of sequences and sequence data from amplicon id 
        
        Remember amplicon id is just one amplicon that uses the primer set. We want to return all
        amplicons that use the primer set of the given amplicon.
        """

        k = super().lookup_param("k")
        seqs = self.dbc.execute("""
                SELECT Genomes.name, Contigs.name, Amplicons.start1 as start, Amplicons.start2 + {k} as end, Amplicons.strand, SUBSTR(Sequences.seq, Amplicons.start1 + 1, Amplicons.start2 + {k} - Amplicons.start1) 
                    
                FROM Amplicons 
                        
                INNER JOIN Genomes 
                    ON Amplicons.genome = Genomes.id 
                
                INNER JOIN Contigs 
                    ON Amplicons.contig = Contigs.id 
                    
                INNER JOIN Sequences 
                    ON Amplicons.contig = Sequences.id 
                            
                WHERE Amplicons.kmer1 = (SELECT kmer1 FROM Amplicons WHERE rowid = '{rowid}' LIMIT 1) AND
                    Amplicons.kmer2 = (SELECT kmer2 FROM Amplicons WHERE rowid = '{rowid}' LIMIT 1)
                        
                """.format(k=k, rowid=amp_id)).fetchall()
 
        return seqs


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
            self.locations.append((genome, contig, start+indx, strand, kmer_id))
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
       

def lazy_imap(processes, function, iterable):
    """ 
     Analog of pool.imap that behaves lazily to save memory 
    
    This exploits queue sizes to try to maximize parallelization while minimizing memory.

    The results queue can only hold one completed result per process on the idea that there is no point in having a large queue for the main thread to read from.

    The input queue to the processes only holds processes * 2 items on the idea that there is no point in having a large input queue if the processes are blocked until the results queue is dumped. This helps to ensure that the input queue doesn't fill up in cases of huge iteration.

    The main thread starts the workers and then primes the input queue with 2X processes tasks. Then it enters a loop until all tasks are done that adds a task if there are any remaining and yields a result.


    If the workers move faster than the main thread, they will be mostly blocked waiting for an open slot in the results queue. 

    If the main thread works faster than the workers, it will be mostly blocked waiting for a result to appear.
    """

    input_queue = multiprocessing.Queue()
    results_queue = multiprocessing.Queue(processes)

    # spawn workers
    workers = []
    for indx in range(processes):
        worker = ImapWorker("ImapWorker{}".format(indx), input_queue, results_queue, function)
        worker.start()
        workers.append(worker)

    # ensure the input type is an iterable 
    batch_iter = iter(iterable)

    # keep track of how many batches we are waiting for
    batches_in_progress = 0

    LOG.debug("Priming batch iterator...")
    for _ in range(processes * 2):
        try:
            input_queue.put(next(batch_iter))
            batches_in_progress += 1
        except StopIteration:
            batch_iter = None
            break

    LOG.debug("Yielding Results...")
    while batches_in_progress != 0:
        # check if there are remaining items to be added to the input
        if batch_iter is not None:
            try:
                input_queue.put(next(batch_iter))
                batches_in_progress += 1
            except StopIteration:
                batch_iter = None

        # wait for a result
        yield results_queue.get()
        batches_in_progress -= 1

    for _ in workers:
        input_queue.put("STOP")

    for worker in workers:
        LOG.debug("Waiting for {} to join...".format(worker.name))
        worker.join()

class ImapWorker(multiprocessing.Process):
    """ A worker for the lazy_imap function """
    
    def __init__(self, name, input_queue, results_queue, function):
        super().__init__()

        self.name = name
        self.input_queue = input_queue
        self.results_queue = results_queue
        self.function = function

    def run(self):
        LOG.debug("Starting process {}".format(self.name)) 
        while True:

            # pull an item off the input queue
            itm = self.input_queue.get()

            # check for poison pill
            if itm == "STOP":
                LOG.debug("Returning from {}".format(self.name))
                return
            else:
                # run some function on the item and put the result in the results queue
                result = self.function(itm)
                self.results_queue.put(result)



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

    #counter.count_kmers(fastas)
    #counter.generate_fuzzy_kmers()
    counter.find_conserved_kmers(fuzzy=False)
    #counter.find_potential_amplicons()
    counter.get_amplicon_stats()

def subcommand_main(args):
    main_db = MainDatabase(args.db, args.k, args.threads, args.mismatches, args.stable)  

    main_db.count_kmers(args.fastas, recalculate=args.force)

    if args.fuzzy:
        main_db.generate_fuzzy_kmers(recalculate=args.force)

    main_db.dbc.close()

def subcommand_run(args):
    run_handler = DatabaseRun(args.db, args.output_dir, args.prefix, args.threads) 

    run_handler.process(args.genomes, args.fuzzy, args.amp_min_len, args.amp_max_len, args.amp_frac_genomes, args.force)

    run_handler.dbc.close()


    

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-db", help="path to a new/existing main database", required=True)
    parser.add_argument("-force", help="force overwriting existing paths", action="store_true")
    parser.add_argument("-threads", help="number of processes to use [%(default)s]", default=1, type=int)
    parser.add_argument("-fuzzy", help="use fuzzy tables with a particular number of mismatches (default is to not use fuzzy kmers because this is much faster", action="store_true")
    parser.add_argument("-mismatches", help="number of mismatches to allow in fuzzy kmers [%(default)s]", default=1)
    parser.add_argument("-stable", help="number of bp at each end to not allow to be fuzzy [%(default)s]", default=2)

    subparsers = parser.add_subparsers()

    parser_main = subparsers.add_parser("main")
    parser_main.set_defaults(func=subcommand_main)
    parser_main.add_argument("-fastas", help="a bunch of fastas to compare", nargs="+", required=True)
    parser_main.add_argument("-k", help="value of k to use", required=True, type=int)


    parser_run = subparsers.add_parser("run")
    parser_run.set_defaults(func=subcommand_run)
    parser_run.add_argument("-output_dir", help="the directory in which to store the outputi [%(default)s]", default=os.getcwd())
    parser_run.add_argument("-prefix", help="a run specific prefix (this will be the base filename for all files [%(default)s]", default="run")
    parser_run.add_argument("-genomes", help="genome names from the main table to include in this run", nargs='+')
    parser_run.add_argument("-amp_frac_genomes", help="minimum fraction of genomes required to be a good amplicon [%(default)s]", default=1)
    parser_run.add_argument("-amp_min_len", help="heuristic minimum amplicon length to narrow results [%(default)s]", default=100)
    parser_run.add_argument("-amp_max_len", help="maximum amplicon length [%(default)s]", default=400)
    args = parser.parse_args()
    args.func(args)
