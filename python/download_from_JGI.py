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


# TODO Add support for OID in addition to project ID ---Done
        # DONE: but works by passing a (fake) user agent string to the JGI CGI. They
        # don't want bots interfacting with that (as per the "no bots 
        # allowed" message) so I don't want to use that route to access the data. 

# TODO Add support for GenBank and other IMG files. ---Done - GBK

# TODO Raise errors when issues arise rather than returning blanks. Handle these using
        # *specific* try: except statements.

# TODO Naming prefix for download files by id


# these globals are the base paths for the JGI webportals that are needed
global LOGIN, XML, DATA, LOOKUP
LOGIN = 'https://signon.jgi.doe.gov/signon/create'
XML = 'http://genome.jgi.doe.gov/ext-api/downloads/get-directory?organism='     # complete with organism name
DATA = 'http://genome.jgi.doe.gov/'     # complete with path to file from the organism XML
LOOKUP = 'http://genome.jgi.doe.gov/lookup?'
IMG_LOOKUP = 'https://img.jgi.doe.gov/cgi-bin/w/main.cgi?'

class JGIEntry(object):

    def __init__(self, parent, id, prefix="", convert=False):
        self.id = id
        self.parent = parent
        self.prefix = prefix
        self.convert = convert
        #print parent
        self.id_type = "proj"       # temp

        self.name = self._lookup_name()
        if not self.name:
            self.tree = None
            return

        self.tree = self.parse_xml()

    def print_available_files(self, files={},  level=0):
        """ Method to allow recursive printing of files 
        
        My filesystem is a hash but the end file has data that is also stored as a hash
        To get around this right now, I'll look for one of their keys, but really I should do this another way
        Perhaps by using the JGI links to set up a virtual file system, now that would be neat.
        """

        if not files:
            print("Files available for {}".format(self.id), end="\n\n")
            files = self.tree

        for name in sorted(files.keys()):
            if "filename" in files[name]:
                print("   " * level + name)
            else:
                print("   " * level + name + "/")
                self.print_available_files(files=files[name], level=level + 1)
        
    def _lookup_name(self):
        """
        Looks up organism name based on project id

        """
        
        my_params = {'keyName': 'jgiProjectId', 'keyValue': self.id}
        name_html = self.parent.session.get(LOOKUP, params=my_params)    # use the parent session to maintain cookies
        print((self.id, name_html.url))
        
        
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
            print("Name not found for id: {}".format(str(self.id)))

            return ""
            # instead of return, might be nice to remove the reference in the parent list

    def lookup_taxon_oid(self):
        """ Coming soon! """
        raise NotImplementedError

    def parse_xml(self):
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

    def download_data(self, download_regex, resume=False, local_dir=os.getcwd(), make_dir=False):
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
        # if there is no tree skip 
        if not self.tree:
            print("No files (or access denied) for {}. Skipping Download.".format(self.id))
            return

        if make_dir:
            local_dir = local_dir + "/" + self.id
            if not os.path.isdir(local_dir):
                try:
                    os.mkdir(local_dir)
                except IOError as e:

                    print(e)
                    sys.exit("Could not make directory")

        
        to_download = self._get_download_tuple(download_regex)

        for url, local_name in to_download:
            if resume:
                if os.path.isfile(local_name):
                    print("File: {} found locally. Skipping download.".format(local_name))
                    continue
            print("Downloading {} to\n    {}".format(url, local_dir + "/" + local_name))
            self._download_path(url, local_dir + "/" + local_name)

    def _get_download_tuple(self, download_regex, concat="."):
        """
        Given download_regex = [(folder1, regex1), (folder2, regex2)]
        this method returns = [(URL1, local_path1), (URL2), local_path2)]

        Creates local path based on prefix + suffix where:
            prefix = one of {self.prefix -> self.name -> self.id} in that order
            suffix = one of {captured group 1 -> full JGI filename} in that order

        Concatenates prefix and suffix together with the string specified by concat

            """
        
        download_files = []
        # build an array of urls to download while making sure everything checks out as expected
        for path, regex in download_regex:
            directory = {}
            parent_dir = "root"
            for folder_re in path.split("/"):
                if not folder_re:
                    continue
                if not directory:
                    directory = self.tree
                for file_name in directory:
                    #print(file_name)
                    match = re.search(folder_re, file_name)
                    if match and "filename" not in directory[file_name]:     # match & not file
                        directory = directory[file_name]
                        parent_dir = file_name
                        break
                else:
                    print("Folder '{}' not found in '{}' for {}".format(folder_re, parent_dir, self.id))
                    directory = ''

            if directory:

                # get files that the regex matches
                new_files = []
                for file_name in directory.keys():
                    match = re.search(regex, file_name)
                    if match:
                        url = directory[file_name]["url"]
                        # this will be 2 step naming - prefix = org name or user supplied
                        # - suffix = capture group or full filename
                        # try for prefix -> organism name -> id
                        if self.prefix:
                            prefix = self.prefix
                        elif self.name:
                            prefix = self.name
                        else:
                            prefix = self.id
                            
                        # make the suffix
                        # option to name it from a captured group
                        if len(match.groups()) >= 1:
                            suffix = match.group(1)
                        else:
                            suffix = file_name

                        # merge the prefix and suffix with "_" to give final filename
                        new_files.append((url, prefix + concat + suffix))
                
                # see if the regex was successful and add files it found to files the other regexes found
                if new_files:
                    download_files += new_files
                else:
                    print("No files found for regex: '{}', id: '{}.'".format(regex, self.id), file=sys.stderr)
                    print("Skipping downloading...", file=sys.stderr, end="\n\n")
                    #return []

            else:
                #return []
                pass
        
        # make sure all local names are unique
        # NOTE: This will still overwrite already existing local names, may want a "force" option
        # to enable default to be to not overwrite files
        names = [name for url, name in download_files]
        if len(names) > len(set(names)):
            print("Attempted to download multiple files using the same local path.", file=sys.stderr)
            for url, name in download_files:
                print("  {}\t{}".format(name, url), file=sys.stderr)
            print("Skipping downloading...", file=sys.stderr, end="\n\n")
            return []

        return download_files

    def _download_path(self, url, local_path):
        url = DATA + url
        # print(url)

        #download the file in chunks
        request = self.parent.session.get(url, stream=True)
        with open(local_path, 'wb') as fh:
            for chunk in request.iter_content(chunk_size=1024):
                if chunk:  # filter out keep-alive new chunks
                    fh.write(chunk)
                    fh.flush()


class IMGEntry(object):
    """ 
    A class to get data from the IMG web portal.

    Note: This class is different than a JGI entry because this one uses taxon ids
    to interact with the CGI on IMG's web portal to download a different set of data.

    This class will effectively help to clean up the JGI entry class by splitting the 
    work.
    """

    DATA_TYPES = ['gbk']


    def __init__(self, parent, taxon_oid, prefix="", convert=False):
        self.parent = parent
        self.taxon_oid = taxon_oid
        self.prefix = prefix
        self.convert = convert

        # header to display to the website
        self.header = {'User-Agent': 'Mozilla/5.0 (X11; Linux x86_64; rv:24.0) Gecko/20140924 Firefox/24.0 Iceweasel/24.8.1'}

        #self.name = self._lookup_name()

    def lookup_proj_id(self):
        """ 
        Returns: the JGI project id as a string.
        Excepts: ValueError
        """
            
        my_params = {'section': 'TaxonDetail', 'page': 'taxonDetail', 'taxon_oid': self.taxon_oid}
        img_html = self.parent.session.get(IMG_LOOKUP, params=my_params, headers=self.header)
        
        # match the keyValue param in the link with digits to allow the keyValue param to
        # be anywhere within the link 
        # this will need to be changed if JGI ever adds letters to their ids
        m = re.search("\"genome-btn download-btn.*href=['\"].*keyValue=(\d*).*['\"]", img_html.text)
        if m:
            print(m.group(1))
            return m.group(1)
        else:
            #print("Couldn't find URL for IMG id: {}".format(self.taxon_oid), file=sys.stderr)
            raise ValueError("Download link for found for IMG Id: {}".format(self.taxon_oid))

    def print_available_files(self):
        """ 
        Prints available files form IMG, if convert enabled, also prints for JGI
        Currently does no checking, just tells the file types that are implemented 
        """

        print("Files available from IMG for taxon oid {}:".format(self.taxon_oid))
        print("\t" + "\n\t".join(self.DATA_TYPES))

        if self.convert:
            jgi = JGIEntry(id=self.lookup_proj_id(), parent=self.parent)
            jgi.print_available_files()
    


    def download_data(self, download_type, resume=False, force=False):
        if download_type in self.DATA_TYPES:
            if download_type == 'gbk':
                gbk_file = self.download_gbk()
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
        
        payload = {     'taxon_oid': self.taxon_oid,
                        'scaffold_oid': 'all',
                        'format': 'gbk',
                        '_section_TaxonDetail_processArtemisFile': 'Go'
                    }
        
        
        # NOTE: this can take a while for genomes with many scaffolds
        # ALSO NOTE: Don't use this! For genomes with > 1000 scaffolds,
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
            print(request.text)
            raise ValueError("PID not found.")
     
        # This is the payload required to download the file. PID obtained from request
        payload = {     'pid': pid,
                        '_section_TaxonDetail_downloadArtemisFile_noHeader': 'Download File',
                        'type': 'gbk'
                    }

        # dowload the .gbk file
        r = self.parent.session.post(IMG_LOOKUP, data=payload, headers=self.header)
    
        # TODO error checking here to make sure gbk_file is actually a gbk file

        return r.text


class JGIInterface(object):
    """
    Class for managing the connection with JGI, responsible for logging in, establishing, and maintaining cookies
    throughout the session
    """

    def __init__(self, login_file):
        self.entries = []
        self.session = self._login(login_file)

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
            print("Login file not specified")
            usage()
            sys.exit(2)

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
            print("Login Successful!", file=sys.stderr)
            return session

    def make_entries(self, ids, id_type, prefix="", convert=False):
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
                entry = JGIEntry(parent=self, id=id, prefix=prefix, convert=convert)

            elif id_type == "img_oid":
                entry = IMGEntry(parent=self, taxon_oid=id, prefix=prefix, convert=convert)
            else:
                raise AssertionError("Wrong id type, should have never gotten here.")

            self.entries.append(entry)

    @staticmethod
    def _read_file_or_split_arg(arg, delim="|"):
        """ 
        Either reads the file from a filename or splits the arg
        Returns: A list of items
        Excepts: None
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


    def download(self, url):
        #download the file in chunks

        local_path = url.split("/")[-1]
        print("Downloading {} to \n{}".format(url, local_path))
        request = self.session.get(url, stream=True)
        with open(local_path, 'wb') as fh:
            for chunk in request.iter_content(chunk_size=1024):
                if chunk:  # filter out keep-alive new chunks
                    fh.write(chunk)
                    fh.flush()


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

def test_download_gbk(interface, img_oid):
    genome = IMGEntry(parent=interface, taxon_oid=img_oid)

    with open("test.gbk", 'w') as OUT:
        OUT.write(genome.download_genbank())


def main(args):
    """
    Main function to manipulate the classes to download data.
    """
    
    interface = JGIInterface(args.login)
    interface.make_entries(args.ids, args.id_type, prefix=args.names, convert=args.convert)

    ##### TEMP
    #test_download_gbk(interface, args.ids)
    #sys.exit()

    # if raw download from path
    if args.download:
        interface.download(args.download)
        sys.exit()

    if args.get:
        for entry in interface.entries:
            entry.download_data(args.get, args.r)
    else:
        for entry in interface.entries:
            entry.print_available_files()
       # interface.entries[0].print_available_files()


if __name__ == "__main__":

    os.nice(20)     # give this program lowest priority to give up CPU to others while downloads are going on

    prefix = os.getcwd()

    # is there a way to remove the whitespace at the beginning of the lines and still have it display correctly?
    # I want this to look like this on the output screen but I want it to be formatted python-style here

    parser = argparse.ArgumentParser(description="""
    Author: Hunter Cameron
    Github: https://github.com/hunter-cameron/Bioinformatics/tree/master/python
    
    Description:
        Downloads files from the Joint Genome Institute (JGI) based on JGI
        Project Id. 

        Files to be downloaded are specified using Python-style regular
        expressions. 

        Don't supply option '--get' to print a list of file available for 
        download.

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
    parser.add_argument("--get", "-g", help="""
folder and regex of data to download separated by '|'. More than 
one can be included as positional arguments separated by a space. 
Ex: -g 'folder|regex' 'folder|regex'

    Sample: -g 'IMG Data|.*(?:fna|faa|gbk)$' 
        -- matches any file that ends in fna, faa, 
        or gbk in the IMG Data folder.

Don't supply this option to display files available for download.

Optionally, files can be named using the first capturing group in the 
regex. In the absence of a first capturing group, the 'prefix' for that
id will be used, in absence of the 'prefix', the organism name of JGI
will be used.
        
    """ , type=arg_type_get, nargs='*')
    parser.add_argument("--convert", "-c", help="""
Lookup alt ids to attempt to download the specified data.

Ex. If ids are proj_ids and --get is "gbk", will convert each
proj_id to taxon_oid to get the data.

""", action="store_true")
    parser.add_argument("-names", help="""
file in the same order of 'ids' that specify a prefix for each
id. Alternatively, can be a list on the command line separated
by '|'

    """, type=str)
    parser.add_argument("-d", "--download", help="""
single file (as URL) to download

    """)
    parser.add_argument("--r", "--resume", help="""
resumes file downloads, skips files that already have a local path

    """, action="store_true")
    args = parser.parse_args()

    main(args)

