

import argparse
import sys
import os
import pandas
import getpass
import shutil
import tarfile
import logging
import subprocess
from Bio import SeqIO

from mypyli import utilities
from mypyli.jgi_interface import JGIInterface, JGIOrganism
logging.basicConfig(level=logging.DEBUG)
LOG = logging.getLogger(__name__)

class MissingDataError(ValueError):
    """ Exception to raise when trying to get the path for data that isn't present. """
    pass

class JGIDownloadError(Exception):
    """ Blanket error to use when JGI download fails for any reason"""
    pass


class DataReport(object):
    """ Custom data type to store info about data"""

    def __init__(self):
        self.found_ids = []
        self.missing_ids = []
        self.extra_ids = []
        self.invalid_files = []
        self.duplicates = []

    def __str__(self):
        return "\n".join([  "Data Report:",
                            "    Found ids: " + str(self.found_ids),
                            "    Missing ids: " + str(self.missing_ids),
                            "    Extra ids: " + str(self.extra_ids),
                            "    Invalid files: " + str(self.invalid_files),
                            "    Duplicates: " + str(self.duplicates)
                            ])

    def merge(cls, rep1, rep2):
        """ Merges two reports and returns the merged version. """
        merged = cls()
        # I could do this merge in a for loop by accessing the internal dict
        # however, I think its frowned upon to mess with the internals
        merged.found_ids = rep1.found_ids + rep2.found_ids   
        merged.missing_ids = rep1.missing_ids + rep2.missing_ids
        merged.extra_ids = rep1.extra_ids + rep2.extra_ids


class DataType(object):
    """ DataType with associated paths and files """

    def __init__(self, name):

        self.name = name
        self.files = {}

    def __str__(self):
        return "DataType {}".format(self.name)

    def update(self):
        """ Updates the files dict for whether or not each file exists """
        for path in self.files:
            self.files[path] = os.path.exists(path)

    def add_data_file(self, path):
        """ 
        Adds a file to the datatype
        
        Prefix is a full path to the file
        """
        
        if path in self.files:
            LOG.warning("Refusing to add duplicate data path '{}' to {}".format(path, str(self)))
            return
        else:
            self.files[path] = False

    def get_missing(self):
        """ Returns a list of missing files """
        return [file for file in self.files if not self.files[file]]

    def get_data_paths(self):
        """ Returns a list of data paths """
        return list(self.files.keys())

    @property
    def present(self):
        """ Property that is True if all the data exists and false otherwise """
        self.update()

        for type_exists_bool in self.files.values():
            if not type_exists_bool:
                return False
        else:
            return True




class Isolate(object):
    """ Representation of an isolate complete with methods to check for existing data and download more.
    
    A major advantage (problem?) of this class is that is reports all failed data download attempts as warnings.
    This is great for interactive use but if you are trying to collect failed downloads using a 
    program, it may be a little more difficult. I could add an option to suppress errors and the
    warning would be logged if errors are suppressed, otherwise the error raised?
    """


    
    DTYPES = ["bundle", "gbk", "genome", "blast_db", "genes", "ko", "cog", "pfam", "tigrfam", "interpro"]
    
    def __init__(self, taxon_oid, name=None):
        """
        Taxon_oid is required for online lookups.

        Name is what you want to name member files as. Defaults to taxon_oid.
        """
        self.taxon_oid = taxon_oid

        # the name attrib will be what datafiles are named
        if name:
            self.name = name
        else:
            self.name = taxon_oid


        # data fields 
        self.datatypes = {}
        self.metadata = {}

        # additional data fields for/from the database
        self.database_data = {}

        # stores a JGIOrganism to manage online data access
        self.organism = None

    def __str__(self):
        return "Isolate: {}".format(self.name)

    def get_data_paths(self, dtype):
        """ Returns a list of the filepath(s) for a particular datatype"""
        return self.datatypes[dtype].get_data_paths()

    def make_organism(self, interface):
        """ Creates a JGIOrganism object from the supplied interface """ 

        self.organism = JGIOrganism(interface, taxon_oid=self.taxon_oid, prefix=self.taxon_oid)

    # methods for looking for data files
    def add_datatype(self, dtype, pre_suf):
        """ Adds a DataType object to the datatypes attribute
        
        datatype - a string name for the datatype
        
        pre_suf - a tuple or list of tuples with a prefix (full directory path) 
        and suffix (extension) for each data file to put in the datatype
        """
       
        try:
            datatype = self.datatypes[dtype]
        except KeyError:
            datatype = DataType(dtype)

        # convert pre_suf to a list is it isn't already
        if isinstance(pre_suf, tuple):
            pre_suf = [pre_suf]

        for tup in pre_suf:
            prefix, suffix = tup

            # add a forward slash to the path if necessary
            if not prefix.endswith("/"):
                prefix += "/"

            # get the full path using the isolate name
            data_path = prefix + self.name + suffix
            datatype.add_data_file(data_path)

        self.datatypes[dtype] = datatype

    def get_missing(self):
        """ Returns a dict with missing datatypes as keys and lists of missing files as values """
        missing_data = {}
        for datatype in self.datatypes.values():
            if not datatype.present:
                missing_data[datatype.name] = datatype.get_missing()

        return missing_data

    # methods for adding new data
    def update_metadata(self, overwrite=False):
        """ Looks up metadata from JGI and stores it as an attribute """

        if self.metadata:
            if overwrite:
                LOG.warning("Overwriting metadata for {}".format(str(self)))
            else:
                LOG.warning("Metadata already present for {}. Skipping lookup. Specify overwrite=True to force.".format(str(self)))
                return

        try:
            self.metadata = self.organism.get_metadata()
        except Exception as e:
            LOG.warning("Could not download metadata for {}.\n\t{}".format(str(self), str(e)))

    def update_data(self, dtype=None, overwrite=False):
        """ Tries to get data according to the dtype; if dtype is none, tries to get all missing data """
        if dtype is None:
            if overwrite:
                dtype = self.datatypes.keys()
            else:
                dtype = self.get_missing().keys()
        else:
            if isinstance(dtype, str):
                dtype = [dtype]

            #
            ## make sure the input dtypes are ok
            #
            filtered_dtype = []
            for d in dtype:
                # check if valid datatypes
                if not d in self.datatypes:
                    LOG.warning("Datatype '{}' does not exist for {}. Omitting.".format(d, str(self)))
                    continue
                
                # check if data is already present
                elif self.datatypes[d].present:
                    if not overwrite:
                        LOG.warning("Datatype '{}' is already present for {}. Omitting. Specify 'overwrite=True' to overwrite.".format(d, str(self)))
                        continue

                filtered_dtype.append(d)

            dtype = filtered_dtype

        #
        ## Start getting the data
        #

        # try to get the bundle first
        if "bundle" in dtype:
            try:
                self._download_from_jgi("bundle")
            except (AttributeError, JGIDownloadError) as e:
                LOG.warning("Could not download 'bundle' for {}.\n\t{}".format(str(self), str(e)))

        # data to get from bundle
        if "genome" in dtype:
            try:
                self._extract_from_bundle("genome")
            except (MissingDataError, IOError) as e:
                LOG.warning("Could not extract 'genome' from 'bundle' for {}.\n\t{}".format(str(self), str(e)))

            
        # data to download from JGI
        for d in ["gbk", "ko", "cog", "pfam", "tigrfam", "interpro"]:
            if d in dtype:
                try:
                    self._download_from_jgi(d)
                except (AttributeError, JGIDownloadError) as e:
                    LOG.warning("Could not download '{}' for {}.\n\t{}".format(d, str(self), str(e)))

        # data I need to derrive from other files
        if "genes" in dtype:
            try:
                self._gbk2faa()
            except MissingDataError as e:
                LOG.warning("Cannot process 'genes' for {}.\n\t".format(str(self), str(e)))

        if "blast_db" in dtype:
            try:
                self._make_blast_db()
            except (MissingDataError, RuntimeError) as e:
                LOG.warning("Cannot 'makeblastdb' for {}.\n\t{}".format(str(self), str(e)))

    def _download_from_jgi(self, dtype):
        """ Downloads a datatype for a given organism """
        # download what we can from JGI, should be try-excepted
        # this doesn't try to change the name of the resultant file
        # need to delete .tar.gz after done

        LOG.debug("Trying to download {} from jgi for {}...".format(dtype, str(self)))
        if dtype == "bundle":
            try:
                self.organism.download_data(("IMG Data", ".*.(tar.gz)$"))
            except Exception as e:
                raise JGIDownloadError("Download failed with error message:\n\t{}".format(str(e)))

            tar = tarfile.open(self.taxon_oid + ".tar.gz", 'r:gz')   

            # we need the directory before the bundle dir to extract into
            # this command joins the path with the parent dir (..) and then resolves an abs path
            extract_dir = os.path.abspath(os.path.join(self.get_data_paths("bundle")[0], os.path.pardir))
            tar.extractall(path=extract_dir)
            os.remove(self.taxon_oid + ".tar.gz")

        # download from IMG
        if dtype in ['gbk', 'ko', 'cog', 'pfam', 'tigrfam', 'interpro']:
            try:
                self.organism.download_data(dtype)
            except Exception as e:
                raise JGIDownloadError("Download failed with error message:\n\t{}".format(str(e)))

            # move the downloaded file to the proper location
            shutil.move(self.taxon_oid + "." + self.organism.IMG_DATA_SUF[dtype],
                        self.get_data_paths(dtype)[0])

    def _extract_from_bundle(self, dtype):
        """ 
        Tries to extract a data type from the bundle. 

        Excepts FileNotFoundError if bundle doesn't exist or doesn't have the file 
        """
 
        # this is a hash that maps dtypes to filenames in the bundle
        # I need to give the user more control over this without digging through the code
        dtype_to_bundle = {
                "genome": self.taxon_oid + ".fna",
                }

        LOG.debug("Trying to extract {} from bundle for {}...".format(dtype, str(self)))
        
            # check if the bundle is available
        if self.datatypes["bundle"].present:
            bundle_path = self.get_data_paths("bundle")[0]
        else:
            raise MissingDataError("Cannot extract data from bundle; bundle has not been downloaded.")

        # look for the file
        for file in os.listdir(bundle_path):
            if file == dtype_to_bundle[dtype]:
                src = bundle_path + "/" + file
                dest = self.get_data_paths(dtype)[0]
                LOG.info("Copying {} to\t{}".format(src, dest))
                try:
                    shutil.copy(src, dest)
                    break
                except Exception as e:
                    raise IOError("Copying {} to {} failed. {}".format(src, dest, str(e))) 
        else:
            raise MissingDataError("Filename {} not found in bundle.".format(dtype_to_bundle[dtype]))

    def _make_blast_db(self):
        """ Makes a blast database using the makeblastdb system command"""

        LOG.debug("Trying to makeblastdb for {}...".format(str(self)))

        # check if the genome is available
        if self.datatypes["genome"].present:
            fasta = self.get_data_paths("genome")[0]
        else:
            raise MissingDataError("'genome' data required for making blast database.")

        # get the prefix
        prefix = self.get_data_paths("blast_db")[0].rsplit(".", 1)[0]
        
        status = subprocess.call(["makeblastdb", 
                        "-dbtype", "nucl", 
                        "-in", fasta,
                        "-out", prefix])
        if status:
            raise RuntimeError("makeblastdb command failed; make sure it is on $PATH")
        else:
            LOG.info("BLAST database successfully created!")

    def _gbk2faa(self):

        LOG.debug("Trying to extract genes for {}...".format(str(self)))

        # check if the gbk is available
        if self.datatypes["gbk"].present:
            gbk = self.get_data_paths("gbk")[0]
        else:
            raise MissingDataError("'genome' data required for making blast database.")

        # get the path for the output
        genes = self.get_data_paths("genes")[0]

        utilities.gbk2faa(gbk, genes)

        
    @property
    def organism(self):
        if self._organism is None:
            raise AttributeError("No JGIOrganism has been created for: {}".format(str(self)))
        else:
            return self._organism
    @organism.setter
    def organism(self, value):
        self._organism = value


class IsolateManager(object):
    """ Manages data files for isolates and keeps an up to date database of which files are downloaded. """

    DATA_DIRS = ["bundle", "gbk", "genome", "blast_db", "genes", "ko", "cog", "pfam", "tigrfam", "interpro"]
    DATA_EXT = {
            "bundle": "/",
            "gbk": ".gbk",
            "genome": ".fna",
            "genes": ".genes.faa",
            "ko": ".ko.txt",
            "cog": ".cog.txt",
            "pfam": ".pfam.txt",
            "tigrfam": ".tigrfam.txt",
            "interpro": ".interpro.txt",
            "blast_db": ".fna.nin"      # this is special becuase it is actually just a prefix
            }

    EXT_TO_DTYPE = {
            "/": "bundle",
            "gbk": "gbk",
            "fna": "genome",
            "genes.faa": "genes",
            "ko.txt": "ko",
            "cog.txt": "cog",
            "pfam.txt": "pfam",
            "tigrfam.txt": "tigrfam",
            "interpro.txt": "interpro",
            "fna.nin": "blast_db",
            "fna.nhr": "blast_db",
            "fna.nsq": "blast_db",
            }

    DEFAULT_DATABASE = "isolate_database"

    def __init__(self, dtype_and_ext, database_path="", base_dir=os.getcwd(), mkdir=False):
        """ 
        Argument dtype_and_ext is required and is a list of tuples of (dtype, extension)

        I'm requiring it in init because I want all the isolates to be initialized with the same data
        """
        
        self.base_dir = base_dir
        self.database_path = database_path
        self.dtype_and_ext = dtype_and_ext
        
        # container to hold Isolate objects 
        self.isolates = {}
        
        # option about whether or not to force the creation of data directories
        self.mkdir = mkdir

        # stores a connection with JGI
        self.jgi_interface = None

        # check for default path and read that db if present, otherwise, prepare to write new one
        if not self.database_path:
            self.database_path = base_dir + "/" + self.DEFAULT_DATABASE

        if os.path.isfile(self.database_path):
            self.read_database()
        else:
            self.create_database()


 
    #######################
    ## Isolate based methods
    #
    
    def create_database(self):
        """ Sets up interface to use a new database """
        LOG.info("Creating new database '{}'".format(self.database_path))

    def read_database(self, taxon_oid_field="taxon_oid"):
        """ Reads a database into Isolate Objects """
        LOG.info("Reading database: {}...".format(self.database_path))
        df = pandas.read_csv(self.database_path, sep="\t")
        df = df.where((pandas.notnull(df)), None)
        #df.fillna("NA", inplace=True)

        # check if the taxon_oid field in in the database
        if taxon_oid_field in df.columns:
            # loop through all the entries in the database and make a dict of the values
            for row in df.index:
                database_dict = {}
                for col in df.columns:
                    database_dict[col] = df.loc[row, col]

                isolate = self.add_isolate(str(database_dict[taxon_oid_field]))

                isolate.database_data = database_dict

    def get_missing(self):
        """ Returns a dict of missing files across all isolates """

        all_data_paths = []
        missing_dict = {}
        for isolate in self.isolates.values():
            isolate_missing = isolate.get_missing()
            for dtype in isolate_missing:
                try:
                    missing_dict[dtype].append(isolate.name)
                except KeyError:
                    missing_dict[dtype] = [isolate.name]

        return missing_dict

    def _make_organisms(self):
        for isolate in self.isolates.values():
            isolate.make_organism(self.jgi_interface)

    def update_data(self, dtype=None, overwrite=False):
        for isolate in self.isolates.values():
            isolate.update_data(dtype, overwrite)

    def update_metadata(self, overwrite=False):
        for isolate in self.isolates.values():
            isolate.update_metadata(overwrite)

    def build_database(self, dtype_to_database, database_to_database, metadata_to_database, order=None, replace=False):
        """ 
        Draws on multiple data sources to generate a database.

        Accepts 3 dicts to add fields to the database:
        
        dtype_to_database - allows the user to rename dtypes to whatever they want
        database_to_database - convert keys read in from previous database to field in this database
        metadata_to_database - convery keys from metadata to a database field
        
        Also accepts an order array for column names where the col names should be the set of all the 
        values in the 3 dicts.

        Database is build using dtype first, then database, then metadata. Only empty fields will be filled.
        This allows the metadata to be populated from the original database rather than having
        to look it up each time. 

        To force metadata fields to replace existing fields, specify replace=True.

        returns a pandas DataFrame
        """



        # begin building dict for database
        db_dict = {}
        for index, isolate in self.isolates.items():
            db_dict[index] = {}

            # add data status - will either be True or False
            for k, dk in dtype_to_database.items():
                try:
                    db_dict[index][dk] = isolate.datatypes[k].present
                except KeyError:
                    LOG.warning("{} didn't have '{}' as a listed datatype.".format(str(isolate), k))
                    db_dict[index][dk] = False

       
            # add data from a previous database -- these should just be manual entry data fields
            for k, dk in database_to_database.items():
                try:
                    db_dict[index][dk] = isolate.database_data[k]
                except KeyError:
                    LOG.warning("{} didn't have '{}' as a database key.".format(str(isolate), k))
                    db_dict[index][dk] = None
 
 
            # add metadata - will either be a string or None
            for k, dk in metadata_to_database.items():
                # skip if the key is already present and not empty and replace isn't specified
                if dk in db_dict[index]:
                    if db_dict[index][dk] is not None and not replace:
                        continue
                
                try:
                    db_dict[index][dk] = isolate.metadata[k]
                except KeyError:
                    LOG.warning("{} didn't have '{}' as a metadata key.".format(str(isolate), k))
                    db_dict[index][dk] = None
        
                        # coerce dict to a DataFrame
        df = pandas.DataFrame.from_dict(db_dict, orient="index")

        # sort the dataframe by the order array
        if order:
            # check for columns in order not in df; this is an error
            for item in order:
                if item not in df.columns:
                    raise ValueError("Column '{}' is present in the order array but not the dataframe.".format(item))
            # check for columns in df not in order; this is a warning
            for item in df.columns.tolist():
                if item not in order:
                    LOG.warning("Column '{}' is present in the df but not in the order array, it will not be in the sorted dataframe.".format(item))


            # return then ordered df
            df = df[order]
            return df

        else:
            return df

    def add_isolate(self, taxon_oid):
        """ 
        Adds an isolate to the dict of isolates this manager manages. 
        
        Returns the isolate for further modification.
        """

        if taxon_oid not in self.isolates:
            isolate = Isolate(taxon_oid)
            
            # add all the data types
            for data_tup in self.dtype_and_ext:
                prefix = self.base_dir + "/" + data_tup[0]
                isolate.add_datatype(data_tup[0], (prefix, data_tup[1]))

            # add the connection to JGI if there is one
            if self.jgi_interface:
                isolate.make_organism(self.jgi_interface)

            self.isolates[taxon_oid] = isolate
            return isolate
        else:
            LOG.warning("Isolate '{}' already in isolate list. Refusing to add duplicate.".format(taxon_oid))

    def connect_to_JGI(self, **kwargs):
        """ 
        Creates a JGIInterface and logs in. 

        Valid kwargs are:
        username -> JGI username
        password -> JGI password
        force_overwrite -> boolean for overwriting existing files
        resume -> boolean for skipping existing files
        newest_only -> boolean for downloading the newest file only when name conflicts occur

        """

        # set kwarg defaults
        kwarg_values = {
                "username": "",
                "password": "",
                "force_overwrite": False,
                "resume": False,
                "newest_only": False
                }

        for key, value in kwargs.items():
            try:
                kwarg_values[key] = value
            except KeyError:
                raise ValueError("'{}' is not a valid kwarg.".format(key))

        if not kwarg_values["username"]:
            kwarg_values["username"] = input("JGI username (email): ")

        if not kwarg_values ["password"]:
            kwarg_values["password"] = getpass.getpass("JGI Password for {}: ".format(kwarg_values["username"]))

        self.jgi_interface = JGIInterface(username=kwarg_values["username"],
                                    password=kwarg_values["password"],
                                    force_overwrite=kwarg_values["force_overwrite"],
                                    resume=kwarg_values["resume"],
                                    newest_only=kwarg_values["newest_only"]
                                    )

        self._make_organisms()

    def check_for_new_isolates(self, projects):
        """ 
        Checks for new isolates in a list of specfied projects. 

        Returns a dict of isolates not in database taxon_oid -> Isolate to allow easy
        lookup of metadata to determine if the isolate should be added.
        """
        # get all the taxon_oids from all the projects
        taxon_oids = set()
        for project in projects:
            tids = self.jgi_interface.get_taxon_oids_for_proposal(project)
            taxon_oids.update(tids)

        # check for new ones
        new_tids = {}
        for taxon_oid in taxon_oids:
            if taxon_oid not in self.isolates:
                isolate = Isolate(taxon_oid, name=taxon_oid)
                isolate.make_organism(self.jgi_interface)
                new_tids[taxon_oid] = isolate 

        return new_tids

    def make_all_blast_db(self):
        """ Concatenates all genomic fastas and makes a blast database out of that. """


        #
        ## Concatenate genomic fastas and prepend name to each header
        #
        LOG.info("Concatenating all files...")
        tmp_fasta_name = "all_isolates.fna"
        with open(tmp_fasta_name, "w") as OUT:
            for isolate in self.isolates.values():
                if isolate.datatypes["genome"].present:
                    # read the fasta file
                    with open(isolate.get_data_paths("genome")[0], "r") as IN:
                        LOG.debug("    Adding file from {}...".format(str(isolate)))
                        for record in SeqIO.parse(IN, "fasta"):
                            # change the header
                            record.id = isolate.name + "_" + record.id

                            # set description to id to avoid double printing
                            record.description = record.id

                            SeqIO.write(record, OUT, "fasta")


        #
        ## Make a mock isolate to access the mkblastdb method
        #
        all = Isolate(None, "all_isolates")
        all.add_datatype("genome", (self.base_dir + "/", ".fna"))

        all.add_datatype("blast_db", [tup for tup in self.dtype_and_ext if tup[0] == "blast_db"])
        all.update_data("blast_db", overwrite=True)

        # remove temporary fasta
        os.remove(tmp_fasta_name)
    
    @staticmethod
    def help():
        """ Displays the help message. """
 
        print("""
Isolate Download Manager initialized as 'manager'.


Other notable variables:
    dtype_to_database -> dict mapping datatypes to field names for tabular output
    metadata_to_database -> dict mapping metadata keys (on IMG organism summary) to field names
    database_to_database -> dict mapping fields in a previous database to field names in this db
    field_order -> list that contains the order of all the fields from the 3 dicts above
    projects -> list of the projects the isolates are currently from

You might want to:

    - Add an isolate -> manager.add_isolate(taxon_oid)

    - See missing files -> manager.get_missing()

    - Connect to JGI -> manager.connect_to_JGI()

    - Get missing data for all isolates -> manager.update_data()
        - Get a specific data type -> manager.update_data("gbk")

    - Get metadata for all isolates -> manager.update_metadata()

    - Get a tabular representation of the isolates -> 
            df = manager.build_database(dtype_to_database, database_to_database, metadata_to_database, field_order)

    - Print the tabular representation -> write_database(df, filename)

    - Check for new isolates in existing projects (doesn't add them to database) ->
            manager.check_for_new_isolates(projects)


    - Display this message again -> manager.help()
""")


def write_database(df, filename):
    """ Writes a pandas DataFrame to the specified filename. """
    df.to_csv(path_or_buf=filename, sep="\t", na_rep='NA', header=True, index=True, index_label="taxon_oid")



        
class UserInterface(object):
    """ Class for a user interface to this program. This is kept for educational purposes only. The method for using this script now is to run it in Ipython. """
    
    STATES = {  "main": 
                    [ 
                        "Welcome to the Isolate Manager",
                        "==============================",
                        "0. Exit", 
                        "1. Update Database",
                        "2. Write Database",
                        "3. Add Isolate"
                    ],

                "AddIsolate":
                    [
                        "Add Isolate from:",
                        "=================",
                        "0. Back",
                        "1. JGI Taxon OID",
                    ],

                "GetData":
                    [
                        "Get Data",
                        "========",
                        "0. Main Menu",
                        "1. Get all", 
                        "2. Get local"
                    ]
            
                    }
                        
                    


    
    RESPONSES = { "main": 
                    [   
                        "Action:Exit",
                        "Action:UpdateDB",
                        "Action:WriteDB",
                        "State:AddIsolate"
                    ],

                "AddIsolate":
                    [
                        "State:main",
                        "Action:AddTaxonId",
                    ],

                "GetData":
                    [
                        "State:main",
                        "Action:GetDataAll",
                        "Action:GetDataLocal"
                    ]

                }

    
    def __init__(self, manager):
        self.manager = manager
        self.running = True
        self.current_state = "main"    


    def display(self):
        """ Main loop of the user interface """
        while self.running:
            print(end="\n\n")
            for prompt in self.STATES[self.current_state]:
                print(prompt)

            print()
            response = input("Please enter your selection: ")

            self.handle_response(response)

        return self.manager

    def handle_response(self, response):
        # make sure response is an integer
        try:
            response = int(response)
        except ValueError:
            # switch to interactive mode
            if response == "i":
                self.running = False
                return
            print("Error: response must be an integer")
            return

        # make sure response is in the list of keys
        try:
            action = self.RESPONSES[self.current_state][response]
        except IndexError:
            print("Error: response must be in the set: " + 
                   " ".join([str(i) for i in range(len(self.RESPONSES[self.current_state]))]))
            return

        # handle the response
        res_type, command = action.split(":")
        
        if res_type == "State":
            self.current_state = command

        elif res_type == "Action":
            # main commands
            if command == "Exit":
                print("Goodbye.")
                sys.exit()
            elif command == "UpdateDB":
                self.update_db()
            elif command == "WriteDB":
                self.write_db()

            # add isolate commands
            elif command == "AddTaxonId":
                self.add_isolate_taxon_oid()
            elif command == "AddProjId":
                self.add_isolate_proj_id()
            elif command == "AddUniqueId":
                self.add_isolate_unique_id()

            # get data commands
            elif command == "GetDataAll":
                self.get_data("all")
            elif command == "GetDataLocal":
                self.get_data("local")

    def update_db(self):
        print()
        self.manager.update_database()
        print()
        print("Database updated.")
        self.current_state = "GetData"

    def write_db(self):
        path = input("Please type the path to write the database (default=" + self.manager.database_path + ").\n")

        if not path:
            path = self.manager.database_path

        self.manager.write_database(path)

        print("Database written to: {}".format(path))


    def add_isolate_taxon_oid(self):
        isolate = input("Type the JGI taxonomy id: ")

        if len(isolate) != 10:
            print("JGI taxonomy ids have 10 digits. Refusing to accept this entry because it does not have 10 characters.")
            return
        else:
            self.manager.add_isolate(isolate)
            print("Isolate {} successfully added.".format(isolate))


    def get_data(self, style):
        if style == "all":

            self.current_state == "main"
            username = input("JGI username (email): ")
            password = getpass.getpass("JGI Password: ")
            jgi_interface = JGIInterface(username=username, password=password, newest_only=True)
            self.manager.get_data(jgi_interface)
        elif style == "local":
            self.manager.get_data()

def get_ids_from_file(ids_f):
    ids = []
    with open(ids_f, 'r') as IN:
        for line in IN:
            ids.append(line.rstrip())

    return ids

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Provides an interactive interface for managing/downloading a JGI/IMG based database. The ideal way to use this is by invoking this script using ipython.")
    parser.add_argument("-db", "-database", help="a database mapping Taxon id to types of data on file. will be created if blank", default="")
    parser.add_argument("-dir", "-base_dir", help="the base isolates directory", default=os.getcwd());
    parser.add_argument("-mkdir", help="make data directories if they don't exist", action="store_true")
    parser.add_argument("-ids", help="list of ids (or a file with a list) to check, this is optional if there is an established database", nargs="*")
    args = parser.parse_args()

    
    args.dir = os.path.abspath(args.dir)

    # make a list of datatype and file extensions
    # this can be tailored to look for less data
    dtype_and_ext = [
            ("bundle", "/"),
            ("gbk", ".gbk"),
            ("genome", ".fna"),
            ("genes", ".genes.faa"),
            ("ko", ".ko.txt"),
            ("cog", ".cog.txt"),
            ("pfam", ".pfam.txt"),
            ("tigrfam", ".tigrfam.txt"),
            ("interpro", ".interpro.txt"),
            ("blast_db", ".fna.nin"),
            ("blast_db", ".fna.nhr"),
            ("blast_db", ".fna.nsq")
            ]

    # options to use when making a database
    field_order = [ "freezer_id",
                    "organism_name",
                    "project",
                    "sequencing_center",
                    "type",
                    "status",
                    "gram_staining",
                    "lineage",
                    "bundle",
                    "gbk",
                    "genome",
                    "genes",
                    "ko",
                    "cog",
                    "pfam",
                    "tigrfam",
                    "interpro",
                    "blast_db"
                    ]

    database_to_database = {
            "freezer_id": "freezer_id",
            "organism_name": "organism_name",
            "sequencing_center": "sequencing_center",
            "type": "type",
            "project": "project",
            "gram_staining": "gram_staining",
            "status": "status",
            "lineage": "lineage"
            }

    metadata_to_database = { 
            "Organism Name": "organism_name",
            "Sequencing Center": "sequencing_center",
            "Culture Type": "type",
            "Study Name (Proposal Name)": "project",
            "Gram Staining": "gram_staining",
            "Lineage": "lineage",
            "Sequencing Status": "status"
            }

    # add all datatypes
    dtype_to_database = {d[0]: d[0] for d in dtype_and_ext} 




    # set up the manager
    manager = IsolateManager(dtype_and_ext, args.db, args.dir, args.mkdir)
    

    field_map = {
            "Organism Name": "organism_name",
            "Sequencing Center:": "sequencing_center",
            "Culture Type": "type",
            "Study Name (Proposal Name)": "project",
            "Gram Staining": "gram",
            }


    # current projects
    projects = {
                'Burkholderia Environmental Isolates',
                'Plant associated metagenomes--Microbial community diversity and host control of community assembly across model and emerging plant ecological genomics systems.',
                'Rhizosphere Grand Challenge Isolate Sequencing'
                }



    # load in any new isolates
    if args.ids is not None:
        for arg in args.ids:
            if os.path.isfile(arg):
                ids = get_ids_from_file(arg)
                for id in ids:
                    manager.add_isolate(id)
            else:
                for id in args.ids:
                    manager.add_isolate(id)

    print("""
Isolate Download Manager initialized as 'manager'.

Other notable variables:
    dtype_to_database -> dict mapping datatypes to field names for tabular output
    metadata_to_database -> dict mapping metadata keys (on IMG organism summary) to field names
    database_to_database -> dict mapping fields in a previous database to field names in this db
    field_order -> list that contains the order of all the fields from the 3 dicts above
    projects -> set of the current projects known to include relevant isolates

You might want to:

    - Add an isolate -> manager.add_isolate(taxon_oid)

    - See missing files -> manager.get_missing()

    - Connect to JGI -> manager.connect_to_JGI()

    - Get missing data for all isolates -> manager.update_data()
        - Get a specific data type -> manager.update_data("gbk")

    - Get metadata for all isolates -> manager.update_metadata()

    - Get a tabular representation of the isolates -> 
            manager.build_database(dtype_to_database, database_to_database, metadata_to_database, field_order)

""")

    #manager.make_organisms(jgi_interface)
    sys.exit()

    ui = UserInterface(manager)
    manager = ui.display()

    #manager.update_database()
