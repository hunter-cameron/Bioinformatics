#!/usr/bin/python

"""
Author: Hunter Cameron
Date: November 4, 2014
Description:

"""

import requests     # not in core python. install using 'pip install --user requests'
import warnings
import re
import os
import sys
import argparse
from argparse import RawTextHelpFormatter
import time

from lxml import etree


# TODO Add some sort of file check to ensure the entire file was downloaded, specifically for the resume option. 

# TODO Genbank file with > 1000 scaffolds -- not sure if this is possible
        # Problem = server-side script checks for scaffold count and demands an email


# these globals are the base paths for the JGI webportals that are needed
global LOGIN, XML, DATA, LOOKUP
LOGIN = 'https://signon.jgi.doe.gov/signon/create'
# I should change this path to use params instead of just completing
XML = 'http://genome.jgi.doe.gov/ext-api/downloads/get-directory?organism='     # complete with organism name
DATA = 'http://genome.jgi.doe.gov/'     # complete with path to file from the organism XML
LOOKUP = 'http://genome.jgi.doe.gov/lookup?'
#IMG_LOOKUP = 'https://img.jgi.doe.gov/cgi-bin/w/main.cgi?'

# This is the url for IMG-ER which might be what I should use because it is the more complete resource according to some folks at JGI.
IMG_LOOKUP = 'https://img.jgi.doe.gov/cgi-bin/er/main.cgi'

# custom exceptions
class DataNotAvailable(Exception):
    """ Exception to be thrown when the requested data cannot be found. """
    pass

class PortalError(Exception):
    """ General exception to be thrown for any type of portal error. """
    pass

class AccessDenied(PortalError):
    """ Exception to be thrown when accessing data fails. """
    pass


class Entry(object):
    """ This is a future parent class for JGIEntry and IMGEntry """
    def __init__(self):
        pass

    def convert_id(self):
        pass


class JGIEntry(object):
    """ Class for looking up data on JGI """
    def __init__(self, parent, id, prefix=""):
        self.id = id
        self.parent = parent
        self.prefix = prefix
        #print parent

        self.taxon_oid = ''      # this will be set in lookup name

        self.name = self._lookup_name()
        if not self.name:
            self.data_tree = None
            return

        self.data_tree = self._parse_xml()

    def _lookup_name(self):
        """
        Looks up organism name based on project id

        At the name lookup step it will be obvious if there are portal errors.

        """
        
        my_params = {'keyName': 'jgiProjectId', 'keyValue': self.id}
        name_html = self.parent.session.get(LOOKUP, params=my_params)    # use the parent session to maintain cookies
        print((self.id, name_html.url))
       


        # go ahead and look up taxon id
        # I know it makes no programming logic sense to do it here but doing it
        # here will prevent another look of this address and be much faster
        taxon_oid_match = re.search('href="https://img.jgi.doe.gov/genome.php?id=(\d)+"', name_html.text)
        if taxon_oid_match:
            self.taxon_oid = taxon_oid_match.group(1)

        # would it be possible to use the url to isolate the name? Rather than the text? I don't think it was becuase I did it this way. I think JGI may not always use the full name in all of their URLs.

        #print(name_html.text)

        # search for href="/(NAME)/ANYTHING.info.html
        match = re.search('href="/(.*)/.*\.info\.html', name_html.text)

        # print(match.group(1))
        # check for a match
        if match:
            self.name = match.group(1)       # get the first matched group
            return match.group(1)
        else:
            raise PortalError("Name lookup failed for {}. URL={}".format(self.id, name_html.url))
            # print("Name not found for id: {}".format(str(self.id)))
            # return ""
            # instead of return, might be nice to remove the reference in the parent list

    # Methods to build data tree
    def _parse_xml(self):
        """
        Parses the XML structure to return a nested hash of folders and files available.

        Structure Expected:

        <organismDownloads>
            <folder1>
                <file1>
                <file2>
            </folder>
            <folder2>
                <file1>
            </folder>
        </organismDownloads>

        Will now handle nested folders.
        """

        link = XML + self.name
        #print "XML:", link
        xml = self.parent.session.get(link).text

        files = {}
        root = etree.fromstring(xml)

        # little sanity check to make sure reading the XML format expected
        assert root.tag == "organismDownloads"

        self._add_children(files, root)

        return files

    def _add_children(self, files, node):
        """
        Recursively adds children nodes to the hash.
        
        This may end up as a class method.
        """

        child_dict = {}
        for child in node:
            if child.tag == "file":
                file_name = child.attrib["filename"]
                child_dict[file_name] = child.attrib

            elif child.tag == "folder":
                folder_name = child.attrib["name"]
                child_dict[folder_name] = {}
                self._add_children(files=child_dict[folder_name], node=child)

            else:
                print("Unknown tag: {}".format(child.tag))

            

        files.update(child_dict)

    # Methods to parse and process URLs to download
    def _get_download_tuple(self, path, file_reg):
        """
        Parses the data tree to get paths that the regex matches.
        Returns: a list of tuples where each entry has the url and the match object from the regex
        """
        
        directory = {}
        parent_dir = "root"

        ### Get the folder the file_regex should be executed in
        for folder_re in path.split("/"):
            # just in case there is a space at the end of the path
            # rstrip may be better here
            if not folder_re:
                continue

            # set the directory if needed (the first loop)
            if not directory:
                directory = self.data_tree

            for file_name in directory:
                #print(file_name)
                match = re.search(folder_re, file_name)
                if match and "filename" not in directory[file_name]:     # match & not file
                    directory = directory[file_name]
                    parent_dir = file_name
                    break
            else:
                #print("Folder '{}' not found in '{}' for {}".format(folder_re, parent_dir, self.id))
                raise DataNotAvailable("Folder '{}' not found in '{}' for {}".format(folder_re, parent_dir,self.id))

        ### Get files that the regex matches
        new_files = []
        for file_name in directory.keys():
            match = re.search(file_reg, file_name)
            if match:
                new_files.append((directory[file_name]["url"], match))
        
        if new_files:
            return new_files
        else:
            #print("No files found for regex: '{}', id: '{}.'".format(file_reg, self.id), file=sys.stderr)
            #print("Skipping downloading...", file=sys.stderr, end="\n\n")
            #return []

            raise DataNotAvailable("No files found for regex: '{}', id: '{}.'".format(file_reg, self.id), file=sys.stderr)

    def _get_full_URL(self, url):
        """ Adds the base path to the URL """
        return DATA + url

    # User Methods
    def print_available_files(self, files={},  level=0):
        """ Method to allow recursive printing of files 
        
        My filesystem is a hash but the end file has data that is also stored as a hash
        To get around this right now, I'll look for one of their keys, but really I should do this another way
        Perhaps by using the JGI links to set up a virtual file system, now that would be neat.
        """

        if not files:
            print("Files available for {}".format(self.id), end="\n\n")
            files = self.data_tree

        for name in sorted(files.keys()):
            if "filename" in files[name]:
                print("   " * level + name)
            else:
                print("   " * level + name + "/")
                self.print_available_files(files=files[name], level=level + 1)

    def lookup_taxon_oid(self):
        """
        This isn't a real lookup because it made more sense to just go ahead and get the taxon oid when parsing out the name but this method is called this for continuity between the two entries (IMG and JGI)
        """

        try:
            return self.taxon_oid
        except AttributeError:
            raise PortalError("Could not find taxon id for {}".format(self.id))

        

    def download_data(self, download_regex):
        """
        Wrapper for _get_download_tuple and _download_path. Allows easier user access.

        make_dir option will make a directory to store all download files from this instance
        
        Example:

            prefix = 090
            XML Structure:
                Raw Data
                    123.myfile.fasta
                    456.yourfile.fasta

            To download all fastas:

                myobj.download_data(download_regex=
                        [("Raw Data", "*.(yourfile.fasta)"), ("Raw Data", "*.(myfile.fasta)")],
                        local_dir=mydir
                        )

            Results:
                mydir/090.myfile.fasta
                mydir/090.yourfile.fasta


            NOTICE: Doing something like:
                myobj.download_data(to_download=[("Raw Data", "*(file).fasta")], mydir)

            will result in an error because there are two files that would be saved to 
            the same name (090.file.fasta)

            See: _get_download_tuple() for more information.

        """
        
        concat = "."
        
        # if there is no tree skip 
        if not self.data_tree:
            print("No files (or access denied) for {}. Skipping Download.".format(self.id))
            return

        # add the files to download
        path, file_reg = download_regex
        to_download = self._get_download_tuple(path, file_reg)

        # make local path prefix
        if self.prefix:
            prefix = self.prefix
        elif self.name:
            prefix = self.name
        else:
            prefix = self.id

        # make the suffix and download the file
        for url, match in to_download:
            
            # make the suffix
            # option to name it from a captured group
            # otherwise name it as the whole match
            if len(match.groups()) >= 1:
                suffix = match.group(1)
            else:
                suffix = match.group()

            # merge the prefix and suffix with "_" to give final filename
            local_name = prefix + concat + suffix
            
            # using the current download method, files could have the same local name
            # which would make then subject to the force option with only one of the
            # multiples actually being written
            self.parent.download_data(self._get_full_URL(url), local_name)


class IMGEntry(object):
    """ 
    A class to get data from the IMG web portal.

    Note: This class is different than a JGI entry because this one uses taxon ids
    to interact with the CGI on IMG's web portal to download a different set of data.
    """

    DATA_TYPES = ['gbk']


    def __init__(self, parent, id, prefix=""):
        self.parent = parent
        self.id = id
        self.prefix = prefix
        
        # header to display to the website
        self.header = {'User-Agent': 'Mozilla/5.0 (X11; Linux x86_64; rv:24.0) Gecko/20140924 Firefox/24.0 Iceweasel/24.8.1'}

        #self.name = self._lookup_name()

    def lookup_proj_id(self):
        """ 
        Returns: the JGI project id as a string.
        Excepts: ValueError
        """
            
        my_params = {'section': 'TaxonDetail', 'page': 'taxonDetail', 'taxon_oid': self.id}
        img_html = self.parent.session.get(IMG_LOOKUP, params=my_params, headers=self.header)
        
        # match the keyValue param in the link with digits to allow the keyValue param to
        # be anywhere within the link 
        # this will need to be changed if JGI ever adds letters to their ids
        m = re.search("\"genome-btn download-btn.*href=['\"].*keyValue=(\d*).*['\"]", img_html.text)
        if m:
            print("Proj id lookup successful. proj_id={}".format(m.group(1)), file=sys.stderr)
            return m.group(1)
        else:
            #print("Couldn't find URL for IMG id: {}".format(self.taxon_oid), file=sys.stderr)
            raise PortalError("Download link not found for IMG Id: {} at URL:\n{}".format(self.id, img_html.url))

    def print_available_files(self):
        """ 
        Prints available files form IMG, if convert enabled, also prints for JGI
        Currently does no checking, just tells the file types that are implemented 
        """

        print("Files available from IMG for taxon oid {}:".format(self.id))
        print("\t" + "\n\t".join(self.DATA_TYPES))
    
    def download_data(self, download_type):
       
        if download_type in self.DATA_TYPES:
            if download_type == 'gbk':
                gbk_file = self.download_genbank()
        else:
            raise ValueError("download_type: {} is not among the acceptable types: {}".format(download_type, str(self.DATA_TYPES)))

    def download_genbank(self):
        """ 
        Returns: a whole genbank file as a string
        Excepts: ValueError

        This function should work but should NOT be used because I have not 
        yet ensured that it will perform properly and have not gotten permission
        from JGI to use this!

        Also, I should probably stream the genbank directly to a file to handle
        the possibility of huge files. This is just a proof-of-concept function
        currently. 
        """
        
        payload = {     'taxon_oid': self.id,
                        'scaffold_oid': 'all',
                        'format': 'gbk',
                        '_section_TaxonDetail_processArtemisFile': 'Go'
                    }
        
        # get the local path
        if self.prefix:
            local_path = self.prefix + ".gbk"
        else:
            local_path = self.taxon_id + ".gbk"

        # check the local path
        if not self.parent.check_local_path(local_path):
            return

        print("\nDownloading genbank file for {}.".format(self.id), file=sys.stderr)

        # NOTE: this can take a while for genomes with many scaffolds
        # ALSO NOTE: Don't use this for genomes with > 1000 scaffolds,
        # JGI has a separate pipeline. I have yet to check if it is the same.

        # this requests the JGI servers generate the file
        request = self.parent.session.post(IMG_LOOKUP, data=payload, headers=self.header)
 
        # get the pid
        # should make this regex more general in case the tags are in a different order
        match = re.search("id='pid' name='pid' value=\'(\d+)\'", request.text)

        try:
            pid = match.group(1)
            print(("pid", pid))
        except:
            #print(request.text)
            
            # attempt to give detailed error mssages
            m = re.search("Please enter your email address since you have selected over 100 entries.", request.text)
            if m:
                raise PortalError("Portal requires email address to download queries with > 100 scaffolds.")

            m = re.search("Please select at least one scaffold", request.text)
            if m:
                raise PortalError("Portal claims no scaffolds were selected (this is a true portal error?).")
            
            
            raise PortalError("PID not found for {}.".format(self.id))
     
        # This is the payload required to download the file. PID obtained from request
        payload = {     'pid': pid,
                        '_section_TaxonDetail_downloadArtemisFile_noHeader': 'Download File',
                        'type': 'gbk'
                    }

        # dowload the .gbk file
        r = self.parent.session.post(IMG_LOOKUP, data=payload, headers=self.header)
    
        # TODO error checking here to make sure gbk_file is actually a gbk file
        
        with open(local_path, 'w') as OUT:
            OUT.write(r.text)

        print("    Downloaded to {}".format(local_path), file=sys.stderr)


class JGIInterface(object):
    """
    Class for managing the connection with JGI, responsible for logging in, establishing, and maintaining cookies
    throughout the session
    """

    def __init__(self, login_file, force=False, resume=False, convert=False):
        self.entries = []
        self.convert_dict = {}

        self.overwritten_paths = []
        self.skipped_paths = []

        self.session = self._login(login_file)

        self.convert = convert
        self.force = force
        self.resume = resume



    @staticmethod
    def _login(login_file=''):
        """
        Parses user credentials from a file and then uses them to set up a session. Returns the session,
        otherwise dies.

        TODO: Need to allow this to accept a username/password rather than a file in case users would rather not make a
        file with their credentials.

        """

        #parse the login info
        if login_file == '':
            raise ValueError("Login file not specified")

        with open(login_file, 'r') as in_handle:
            login_data = {
                'login': in_handle.readline().strip(),
                'password':  in_handle.readline().strip()
            }


        # make a session to store cookie information like a browser would
        session = requests.session()

        # login
        request = session.post(LOGIN, data=login_data)

        # check status of logon
        match = re.search('You have signed in successfully.', request.text)
        if not match:
            print(request.text)       # uncomment to display raw login info. Useful to look through if having problems
            raise AssertionError("Login may not have been successful")
        else:
            print("Login Successful!", file=sys.stderr, end="\n\n")
            return session

    def make_entries(self, ids, id_type, prefix=""):
        """
        Makes a JGIEntry or IMGEntry object for each id in the id_list.
        """

        entries = []

        # get ids
        id_list = self._read_file_or_split_arg(ids)
        
        # get prefixes
        if prefix:
            prefix_list = self._read_file_or_split_arg(prefix)
        else:
            prefix_list = [""] * len(id_list)

        # make all the *entry objects
        for id, prefix in zip(id_list, prefix_list):
            if id_type == "proj_id":
                entry = JGIEntry(parent=self, id=id, prefix=prefix)

            elif id_type == "taxon_oid":
                entry = IMGEntry(parent=self, id=id, prefix=prefix)
            else:
                raise AssertionError("Wrong id type, should have never gotten here.")

            self.entries.append(entry)

    @staticmethod
    def _read_file_or_split_arg(arg, delim="|"):
        """ 
        Either reads the file from a filename or splits the arg
        Returns: A list of items
        """
        items = []
        try:
            with open(arg, 'r') as IN:
                for line in IN:
                    item = line.strip()
                    items.append(item)
        except FileNotFoundError:
            arg = arg.strip()
            items = arg.split(delim)
            
        return items

    def download_data(self, url, local_path):
        """ 
        Downloads a file in chunks.
        Returns: Nothing
        Excepts: AssertionError
        """

        #download the file in chunks

        if self.check_local_path(local_path):

            print("Downloading {} to\n    {}".format(url, local_path), end="\n\n")
            request = self.session.get(url, stream=True)
            with open(local_path, 'wb') as fh:
                for chunk in request.iter_content(chunk_size=1024):
                    if chunk:  # filter out keep-alive new chunks
                        fh.write(chunk)
                        fh.flush()

    def check_local_path(self, local_path):
        """ 
        Checks the local path to ensure no unintentional clobbering.
        
        Note, clobbering can still occur between the time of check and the time of 
        write if there are unknown race conditions.
        
        Nothing in this code should create a race condition.
        """
        if os.path.exists(local_path):
            if self.force:
                print("!!! Force overwriting: {}".format(local_path), file=sys.stderr)
                self.overwritten_paths.append(local_path)
                return True
            else:
                if self.resume:
                    print("Found file: {}. Skipping Download!".format(local_path), file=sys.stderr)
                    self.skipped_paths.append(local_path)
                    return False
                    #raise NotImplemented("This feature will be implemented soon!")
                else:
                    print("Neither force overwrite nor resume are enabled. Refusing to download to {}".format(local_path), file=sys.stderr)
                    raise FileExistsError("Local path must not exist or force/resume must be enabled!")

        else:
            return True





    def convert_entry(self, entry):
        """ 
        Converts an IMGEntry to a JGIEntry and vice versa.
        Returns: Entry object
        Excepts: TypeError
        """

        print("Converting id: {}".format(entry.id), file=sys.stderr)
        if entry in self.convert_dict:
            return self.convert_dict[entry]

        else:
            if type(entry) == JGIEntry:
                conv_entry = IMGEntry(parent=self, id=entry.lookup_taxon_oid(), prefix=entry.prefix)
                self.convert_dict[entry] = conv_entry       # store for future quick conversion
                print("Proj id: {} successfully converted to taxon oid: {}".format(entry.id, conv_entry.id), file=sys.stderr)
                return conv_entry
            elif type(entry) == IMGEntry:
                conv_entry = JGIEntry(parent=self, id=entry.lookup_proj_id(), prefix=entry.prefix)
                self.convert_dict[entry] = conv_entry
                print("Taxon oid: {} successfully converted to proj id: {}".format(entry.id, conv_entry.id), file=sys.stderr)
                return conv_entry
            else:
                raise TypeError("Object: {} is of the wrong type.".format(entry))


def arg_type_get(arg):
    """
    Turns a text string input from the command line to a list of touples like the program requires.
    """

    # see if get is for IMG data, will check for id type later
    if arg in IMGEntry.DATA_TYPES:
        return arg

    try:
        (folder, regex) = arg.split("|", 1)
        return folder, regex
    except:
        raise
        #raise argparse.ArgumentTypeError("Error processing the --get argument.")

def main(args):
    """
    Main function to manipulate the classes to download data.
    
    Maybe should be incorporated into JGIInterface.
    """
    
    interface = JGIInterface(args.login, force=args.force, convert=args.convert, resume=args.resume)
    interface.make_entries(args.ids, args.id_type, prefix=args.names)


    # if download from path
    # perhaps let user specify name for file using the names argument?
    if args.download:
        local_path = args.download.split("/")[-1]
        interface.download_data(args.download, local_path)
        sys.exit()
    
    if args.get:
        failed_log = []
        # loop through entries collecting failed
        for entry in interface.entries:
           
            # loop through regexes collecting failed
            failed_regex = []
            for download_regex in args.get:
                # check if the data types match the entry
                if (type(download_regex) == tuple and type(entry) == JGIEntry) or \
                            (type(download_regex) == str and type(entry) == IMGEntry):
                    download_entry = entry

                else:
                    if args.convert:
                        download_entry = interface.convert_entry(entry)
                    else:
                        raise ValueError("Incorrect --get type for your id. Specify --convert to convert between ids.")

                # try to download data
                try:
                    download_entry.download_data(download_regex)
            
                except (DataNotAvailable, PortalError) as e:
                    print(e, file=sys.stderr)
                    failed_regex.append(download_regex)
            
            # write to the failed log
            failed_log += ["{}: failed data download '{}'".format(entry.id, regex) for regex in failed_regex]

        if failed_log:
            print("\n\nLog of failed downloads:")
            [print("    " + entry) for entry in failed_log]
        else:
            print("\n\nAll files downloaded successfully!")
        if interface.overwritten_paths:
            print("\n\nOverwritten paths:")
            [print("    " + path) for path in interface.overwritten_paths]

        if interface.skipped_paths:
            print("\n\nSkipped paths:")
            [print("    " + path) for path in interface.skipped_paths]
        
    else:
        # this prints available files for all ids, sometimes this may not be wanted
        # right now, I'll just give a message telling how to interrupt
        for entry in interface.entries:
            print("Printing available files for all ids. Press Ctrl+c to interrupt.", file=sys.stderr)
            entry.print_available_files()

            if args.convert:
                conv_entry = interface.convert_entry(entry)
                conv_entry.print_available_files()

if __name__ == "__main__":

    os.nice(20)     # give this program lowest priority to give up CPU to others while downloads are going on

    prefix = os.getcwd()

    # is there a way to remove the whitespace at the beginning of the lines and still have it display correctly?
    # I want this to look like this on the output screen but I want it to be formatted python-style here

    parser = argparse.ArgumentParser(description="""
    Author: Hunter Cameron
    Github: https://github.com/hunter-cameron/Bioinformatics/tree/master/python/
    
    Description:
        Downloads files from the Joint Genome Institute (JGI) based on JGI
        Project id or IMG Taxon Original Id.

        Files to be downloaded are specified using Python-style regular
        expressions. See --get option documentation for more details.

        Don't supply option '--get' to print a list of files available for 
        download. 
        NOTE: this will begin printing the files available for
        every id in the list. If there are many ids, this can take quite some
        time. Interrupt at any time using Ctrl+C.

    """, formatter_class=RawTextHelpFormatter)
    parser.add_argument("--login", "-l", help="""
file with login on one line and password on next

    """, type=str, required=True)
    parser.add_argument("-ids", help="""
file with JGI OIDs, one per line. Alternatively can be supplied 
as a single id or a list of ids separated by '|' at the command prompt
            
    """, type=str)
    parser.add_argument("-id_type", help="""
the type of id supplied

    """, choices=['proj_id', 'taxon_oid'], required=True)
    parser.add_argument("-names", help="""
file in the same order of 'ids' that specify a prefix for each
id. Alternatively, can be a list on the command line separated
by '|'.

    """, type=str)
    parser.add_argument("--get", "-g", help="""
To download data stored on JGI:

Supply the folder and regex of data to download separated by '|'. 
More than one can be included as positional arguments separated 
by a space. Ex: -g 'folder|regex' 'folder|regex'

    Sample: -g 'IMG Data|.*(?:fna|faa|gbk)$' 
        -- matches any file that ends in fna, faa, 
        or gbk in the IMG Data folder.

Don't supply this option to display files available for download.

Optionally, files can be named using the first capturing group in the 
regex. In the absence of a first capturing group, the whole regex match
for that download will be used.

----------------------------------------------------------

To download data stored on IMG, just write the data type.

    Sample: -g 'gbk' 
        -- downloads the GenBank file for an entry
        
    """ , type=arg_type_get, nargs='*')
    parser.add_argument("--convert", "-c", help="""
Lookup alt ids to attempt to download the specified data.

Ex. If ids are proj_ids and --get is "gbk", will convert each
proj_id to taxon_oid to get the data. Otherwise, program will exit
if IMG is attempted to be downloaded from JGI.

When listing files, provide this option to list files from both IMG and 
JGI.

""", action="store_true")
    parser.add_argument("-d", "--download", help="""
single path (as URL) to download

    """)
    parser.add_argument("-f", "--force", help="""
forces download of files, overwrites files that already have a local path

otherwise, program will not overwrite paths.

    """, action="store_true")
    parser.add_argument("-r", "--resume", help="""
skips downloading if local file is found

if this option is not given and a duplicate file is encountered
without the --force option, program will exit

file must have the same name to be considered whole file.

eventually plan to implement size or MD5 checks to ensure
the file is complete before skipping

    """, action="store_true")
    
   
    args = parser.parse_args()

    main(args)

