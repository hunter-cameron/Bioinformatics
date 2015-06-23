

import argparse
import sys
import os
import pandas
import getpass
import shutil
import tarfile

from mypyli import utilities
from mypyli.download_from_JGI import JGIInterface


class MissingDataError(ValueError):
    """ Exception to raise when trying to get the path for data that isn't present. """
    pass

class BundleError(Exception):
    pass


class DataReport(object):
    """ Custom data type to store info about stored data. """

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



class Isolate(object):

    def __init__(self, taxon_oid):
        # JGI information
        self.taxon_oid = name
        self.proj_id = None
        self.organism_name = None

        self.freezer_id = None
        self.taxonomy = None

        # Paths to data
        self.gbk = None
        self.genome = None
        self.genes = None
        self.ko = None
        self.cog = None
        self.pfam = None
        self.tigrfam = None
        self.interpro = None




class IsolateManager(object):
    """ Manages data files for isolates and keeps an up to date database of which files are downloaded. """

    DATA_DIRS = ["bundle", "gbk", "genome", "genes", "ko", "cog", "pfam", "tigrfam", "interpro"]
    DATA_EXT = {
            "bundle": "/",
            "gbk": ".gbk",
            "genome": "fna",
            "genes": ".genes.faa",
            "ko": ".ko.tsv",
            "cog": ".cog.tsv",
            "pfam": ".pfam.tsv",
            "tigrfam": ".tigrfam.tsv",
            "interpro": ".interpro.tsv"
            }

    DEFAULT_DATABASE = "isolate_DB"

    def __init__(self, database_path, base_dir):
        self.database_path = database_path
        self.database = None
        self.isolates = []
        self.base_dir = base_dir
        

        # data information
        self.missing_data = {}

        # check for default path and read that db if present, otherwise, prepare to write new one
        if self.database_path:
            self._read_database()
        elif os.path.isfile(base_dir + "/" + database_path):
            self.database_path = base_dir + "/" + database_path
            self._read_database()
        else:
            self.database_path == base_dir + "/" + database_path


    def _read_database(self):
        """ Reads a database into memory """
        self.database = pandas.read_csv(self.database_path, sep="\t")
        self.database.fillna("NA", inplace=True)

        if "taxon_oid" in self.database.columns:
            self.database.set_index(self.database["taxon_oid"].astype(str), inplace=True)
            self.database.drop("taxon_oid", axis=1, inplace=True)            
            self.isolates = list(self.database.index)
            

    def write_database(self, path):
        self.database.to_csv(path, sep="\t", na_rep="NA", index_label="taxon_oid")


    def add_isolate(self, taxon_id):
        if taxon_id not in self.isolates:
            self.isolates.append(taxon_id)

    def update_database(self):

        # if there is no database, create one
        if self.database is None:
            self.database = pandas.DataFrame(index=self.isolates)

        # add any extra isolates to database if necessary
        for isolate in self.isolates:
            if isolate not in self.database.index:
                new_isolate = pandas.Series([None] * len(self.database.columns), index=self.database.columns, name=isolate)
                self.database = self.database.append(new_isolate)
        
        print(self.database.head())

        missing_data = {}
        extra_data = {}
        potential_duplicates = {}
        non_data = []
        for dir in self.DATA_DIRS:
            # add the dir for the data if is isn't there
            if dir not in self.database.columns:
                self.database[dir] = None

            report = self.check_for_files(dir)
           
            # update the database
            for isolate in report.found_ids:

                #print(isolate in list(self.database.index))
                #print(str(isolate))
                if isolate in self.isolates:
                    self.database.set_value(isolate, dir, "yes")

            #
            ## get some summary info from the update
            #

            # collect missing data types for each isolate
            for isolate in report.missing_ids:
                try:
                    missing_data[isolate].append(dir)
                except KeyError:
                    missing_data[isolate] = [dir]
                
            # collect information about extra isolates (that may should be added?)
            for isolate in report.extra_ids:
                try:
                    extra_data[isolate].append(dir)
                except KeyError:
                    missing_data[isolate] = [dir]

            # report any potential duplicates
            for tup in report.duplicates:
                try:
                    potential_duplicates[dir].append(tup)
                except KeyError:
                    potential_duplicates[dir] = [tup]

            # report non-data files
            for file in report.invalid_files:
                non_data.append(file)



        self.missing_data = missing_data
        
        #
        ## Print a summary 
        #

        print("Database Summary")
        print("================", end="\n\n")

        print("Isolates in database: " + str(len(self.isolates)), end="\n\n")

        print("Isolates with missing data: " + str(len(missing_data)))
        for isolate, datas in missing_data.items():
            print("    " + isolate + ": " + str(datas))
        print()

        print("Isolates not in database with data found: " + str(len(extra_data)))
        for isolate, datas in extra_data.items():
            print("    " + isolate + ": " + str(datas))
        print()

        print("Potentially duplicated files: " + str(len(potential_duplicates)))
        for data_type, tups in potential_duplicates.items():
            print("    " + data_type + ":")
            for tup in tups:
                print("        " + str(tup))
        print()

        print("Non-data files: " + str(len(non_data)))
        for file in non_data:
            print("    " + file)

    def check_for_files(self, dir):
        """ 
        searches a directory for downloaded files and looks for duplicates
        returns a DataReport.
        """
        if dir not in self.DATA_DIRS:
            raise ValueError("Directory must be one of: " + str(DATA_DIRS))


        file_sizes = {}

        report = DataReport()
        for file in os.listdir(self.base_dir + "/" + dir):
            # print(file)
            file_sizes[file] = os.path.getsize(self.base_dir + "/" + dir + "/" + file)
            # All data files should be structured "isolate.datatype.anythingelse"
            isolate, data = file.split(".")[0:2]

            # check if the data is of the right type
            if data != dir:
                report.invalid_files.append(file)
            else:
            
                # check if the isolate was on the list or not
                if isolate in self.isolates:
                    report.found_ids.append(isolate)
                else:
                    report.extra_ids.append(isolate)

        # document any missing ids
        for isolate in self.isolates:
            if isolate not in report.found_ids:
                report.missing_ids.append(isolate)

        # check for duplicates by file size
        # quick check for dups
        # print(file_sizes)
        if len(list(file_sizes.values())) > len(set(file_sizes.values())):
            # longer check to determine which are possible duplicates
            for k1, v1 in file_sizes.items():
                for k2, v2 in file_sizes.items():
                    if v1 == v2:
                        # elim self
                        if k1 != k2:
                            # if converse already entered, skip
                            if (k2, k1) not in report.duplicates:
                                report.duplicates.append((k1, k2))

        return report

    def get_data_path(self, isolate, datatype, check_exists=True):
        """ 
        Gets the path to a specific data file for a specific isolate.

        Set check_exists to false to disable checking if the file exists.
        (useful for generating a path to store data)
        """

        if isolate not in self.isolates or datatype not in self.DATA_DIRS:
            raise ValueError("Isolate or file type does not exist.")

        if check_exists:
            if self.database.get_value(isolate, datatype) == None:
                raise MissingDataError("Datatype {} missing for isolate {}".format(datatype, isolate))
       
        print((isolate, datatype))
        print(self.base_dir + "/" + datatype + "/" + isolate + self.DATA_EXT[datatype]) 
        return self.base_dir + "/" + datatype + "/" + isolate + self.DATA_EXT[datatype] 


    #
    ## Functions to add data
    #

    def get_data(self, jgi_interface=None, missing_data=None):
        """ Adds data as specified by the to_add dict """

        if missing_data is None:
            missing_data = self.missing_data
            
        if jgi_interface:
            
            #
            ## start with remote data
            #

            # make handle for each of the taxon oids
            jgi_interface.make_entries("|".join(missing_data), "taxon_oid", "|".join(list(missing_data.keys())))

            # loop through each isolate
            for entry in jgi_interface.entries:
                # loop through each missing file
                missing = missing_data[entry.id]
                for dtype in missing:
                    
                    # download what we can from JGI, should be try-excepted
                    # this doesn't try to change the name of the resultant file
                    # need to delete .tar.gz after done
                    if dtype == "bundle":
                        centry = jgi_interface.convert_entry(entry)
                        centry.download_data(("IMG Data", ".*.(tar.gz)$"))
                        tar = tarfile.open(entry.id + ".tar.gz", 'r:gz')                        
                        tar.extractall(path=self.base_dir + "/" + "bundle")


                    if dtype in ['gbk', 'ko', 'cog', 'pfam', 'tigrfam', 'interpro']:
                        entry.download_data(dtype)

                        # move the downloaded file to the proper location
                        shutil.move(entry.id + self.DATA_EXT[dtype],
                            self.get_data_path(entry.id, dtype, check_exists=False))

        else:
            #
            ## add local data
            #

            for isolate, missing in missing_data.items():
                for datatype in ["genes"]:
                    if datatype in missing:
                        if datatype == "genes":
                            self._gbk2faa(self.get_data_path(isolate, "gbk"), 
                                    self.get_data_path(isolate, "genes", False))

    def _extract_from_bundle(self, isolate, datatype):
        """ 
        Tries to extract a data type from the bundle. 

        Excepts FileNotFoundError if bundle doesn't exist or doesn't have the file 
        """
        
        datatype_to_bundle = {
                "fasta": 


        bundle = self.get_data_path(isolate, datatype, check_exists=True)
         


    def _gbk2faa(self, genbank, out_path):
        utilities.gbk2faa(genbank, out_path)


def get_ids_from_file(file):
    ids = []
    with open(file, 'r') as IN:
        for line in IN:
            ids.append(line[:-1])

    return ids
        
        
class UserInterface(object):
    
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
                        "2. JGI Project ID",
                        "3. Unique Identifier (data must be downloaded manually)"
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
                        "Action:AddProjId",
                        "Action:AddUniqueId"
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

    def handle_response(self, response):
        # make sure response is an integer
        try:
            response = int(response)
        except ValueError:
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


    def add_isolate_proj_id(self):
        raise NotImplementedError("This function is currently not allowed.");

    def add_isolate_unique_id(self):
        raise NotImplementedError("This function is currently not allowed.");

    def get_data(self, style):
        if style == "all":
            self.current_state == "main"
            username = input("JGI username (email): ")
            password = getpass.getpass("JGI Password: ")
            jgi_interface = JGIInterface("{}\n{}".format(username, password), convert=True)
            self.manager.get_data(jgi_interface)
        elif style == "local":
            self.manager.get_data()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-db", "-database", help="a database mapping Taxon id to types of data on file. will be created if blank")
    parser.add_argument("-dir", "-base_dir", help="the base directory isolates directory", required=True);
    parser.add_argument("-ids", help="list of ids (or a file with a list) to check, this is optional if there is an established database", nargs="*")
    args = parser.parse_args()


    # set up the manager
    manager = IsolateManager(args.db, args.dir)


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

    ui = UserInterface(manager)
    ui.display()

    #manager.update_database()
