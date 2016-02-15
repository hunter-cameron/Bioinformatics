
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

THIS PROGRAM DOES NOT CURRENTLY WORK!!!!


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
