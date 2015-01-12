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
import time

from lxml import etree


# TODO Add support for OID in addition to project ID


# these globals are the base paths for the JGI webportals that are needed
global LOGIN, XML, DATA, LOOKUP
LOGIN = 'https://signon.jgi.doe.gov/signon/create'
XML = 'http://genome.jgi.doe.gov/ext-api/downloads/get-directory?organism='     # complete with organism name
DATA = 'http://genome.jgi.doe.gov/'     # complete with path to file from the organism XML
LOOKUP = 'http://genome.jgi.doe.gov/lookup?'
OID_LOOKUP = 'https://img.jgi.doe.gov/cgi-bin/w/main.cgi?'
#OID_LOOKUP = 'https://img.jgi.doe.gov/cgi-bin/w/main.cgi?section=TaxonDetail&page=taxonDetail&taxon_oid='

def usage():
    print("")
    print("Usage: automate_jgi.py -l <file> -i <file> -o <file> -g 'JGI|raw' (--login, --ids, --out, --get)")
    print("         -l, --login <file>          File with login credentials.(username\\npassword)")
    print("         -i, --ids <file>            File with a list of JGI Taxon Ids")
    print("         -o, --out <file>            Prefix path for output")
    print("         -g, --get <'JGI|raw'>       Folder to download. Currently supports either JGI data or raw data.")
    print("")


class JGIEntry:

    def __init__(self, parent, id, id_type):
        self.id = id
        self.id_type = id_type
        self.alt_id = None      # This will be used to store whichever id is not given
                                # Ex. this is proj_id when taxon_oid is given and vice-versa
                                # This is a first step for allowing files to be named by taxon_oid

        self.downld_prefix = None
        self.parent = parent
        #print parent
        
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
        Looks up organism name based on project Id, would like to add support for IMG taxon oid, but I don't know the keyname.

        Potential method:
            Look up taxon id in IMG using standard link (from website and insert oid)
            Find the download data text. This gives a lookup address identical to what I would have created using
            projid
        """
        
        if self.id_type == 'img_oid':
            # THIS METHOD DOES NOT WORK BECAUSE BOTS DON'T HAVE PERMISSION
            # simulate header to pretend not to be bot
            headers = {'User-Agent': 'Mozilla/5.0 (X11; Linux x86_64; rv:24.0) Gecko/20140924 Firefox/24.0 Iceweasel/24.8.1'}

            
            my_params = {'section': 'TaxonDetail', 'page': 'taxonDetail', 'taxon_oid': self.id}
            img_html = self.parent.session.get(OID_LOOKUP, params=my_params, headers=headers)
            
            m = re.search("\"genome-btn download-btn.*href=['\"](.*)['\"]", img_html.text)
            if m:
                print(m.group(1))
                name_html = self.parent.session.get(m.group(1))
            else:
                print("Couldn't find URL for IMG id: {}".format(self.id), file=sys.stderr)
                return ""

        else:   # id_type = proj
            my_params = {'keyName': 'jgiProjectId', 'keyValue': self.id}
            name_html = self.parent.session.get(LOOKUP, params=my_params)    # use the parent session to maintain cookies
        print((self.id, name_html.url))
        # would it be possible to use the url to isolate the name? Rather than the text It must not be becuase I did it
        # this way. I think JGI may not always use the appropriate name in their URLs.

        #print(name_html.text)
        match = re.search('href="/(.*)/.*\.info\.html', name_html.text)     # search for
        # href="/(NAME)/ANYTHING.info.html

        # print(match.group(1))
        # check for a match
        if match:
            self.name = match.group(1)       # get the first matched group
            return match.group(1)
        else:
            print("Name not found for id: {}".format(str(self.id)))

            return ""
            #instead of return, might be nice to remove the reference in the parent list

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
        """

        link = XML + self.name
        #print "XML:", link
        xml = self.parent.session.get(link).text

        files = {}
        root = etree.fromstring(xml)

        # little sanity check to make sure reading the XML format expected
        assert root.tag == "organismDownloads"

        self._add_children(files, root)

        #for child in root:
        #    if child.tag == "folder":

        #        folder_name = child.attrib["name"]      # need to check if this is correctly assigned
        #        files[folder_name] = {}

        #        for gchild in child:
        #            if gchild.tag == "file":
        #                file_name = gchild.attrib["filename"]
        #                files[folder_name][file_name] = gchild.attrib
        #            else:
        #                print("Unknown grand-child tag: {}".format(gchild.tag))
        #    else:
        #        print("Unknown child tag: {} ".format(child.tag))

        return files

    def _add_children(self, files, node):
        """
        Recursively adds children nodes to the hash
        
        This may end up a class method
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

    def _get_download_tuple(self, download_regex):
        """
        Given a list of files to download in a tuple form (folder, file_regex), this downloads the files and stores
        them using the id as the prefix and and captured groups from the regex as the rest of the name. Returns a
        tuple of (url, local_name)

        Example:

            Genome id = 090
            XML Structure:
                Raw Data
                    123.myfile.fasta
                    456.yourfile.fasta
                    0.file.gbk

            To download all fastas:

                myobj.download_data(to_download=[("Raw Data", "*.(yourfile.fasta)"), ("Raw Data", "*.(myfile.fasta)")],
                save_dir=mydir

            Results:
                mydir/090_myfile.fasta
                mydir/090_yourfile.fasta


            NOTICE: Doing something like:
                myobj.download_data(to_download=[("Raw Data", "*.fasta")], mydir)

            will result in an error because there are two files that could be downloaded.


            This behavior, though unintuitive for people familiar with UNIX wildcards, is intentional to allow the
            user a degree of flexibility in how files are named to make parsing though the downloaded files a bit
            easier. Especially since, currently, JGI bundles files users may be temped to use
            a wild card to collect into a single .tar.gz

            This script is primarily expected to be used to extract a single file from multiple genomes rather than
            many files from a few genomes.

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
                        # option to name it from a captured group
                        if len(match.groups()) >= 1:
                            local_name = self.id + "_" + match.group(1)
                            new_files.append((url, local_name))
                        # otherwise, I'll come up with one. It will be the name + underscore + full filename
                        else:
                            new_files.append((url, self.name + "_" + file_name))

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


class JGIInterface(object):
    """
    Class for managing the connection with JGI, responsible for logging in, establishing, and maintaining cookies
    throughout the session
    """

    def __init__(self, login_file, ids, id_type):

        self.session = self._login(login_file)
        if ids:
            self.genomes = self._make_genomes(ids, id_type)

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
        # print r.text       # uncomment to display raw login info. Useful to look through if having problems
        match = re.search('You have signed in successfully.', request.text)
        if not match:
            raise AssertionError("Login may not have been successful")
        else:
            print("Login Successful!", file=sys.stderr)
            return session

    def _make_genomes(self, ids, id_type):
        """
        Makes a JGIEntry object for each id in the id_list.

        TODO: Accept a list of ids as an option
        """


        genomes = []
        
        # TODO Change this to a try except
        if os.path.isfile(ids):
            #read the jgi_file for proj_ids

            with open(ids, 'r') as in_handle:
                for line in in_handle:
                    line = line.strip()
                    proj_id = JGIEntry(self, line, id_type=id_type)
                    genomes.append(proj_id)

        else:       # parse the ids given directly at the command line
            ids = ids.strip()
            elems = ids.split("|")
            genomes = [JGIEntry(self, genome, id_type=id_type) for genome in elems]


        return genomes

        x = 1
        for genome in self.JGI_genomes:
            print("Downloading Genome: %s" % x)
            x += 1
            genome.download_data("IMG", prefix)

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
    try:
        (folder, regex) = arg.split("|", 1)
        return folder, regex
    except:
        raise argparse.ArgumentTypeError("Error processing the --get argument.")




if __name__ == "__main__":

    os.nice(20)     # give this program lowest priority to give up CPU to others while downloads are going on

    prefix = os.getcwd()

    parser = argparse.ArgumentParser(description="""
    Downloads files from JGI based on Project ID.

    Leave option --get blank to list files available for download

    """)
    parser.add_argument("--login", "-l", help="file with login on one line and password on next", type=str, required=True)
    parser.add_argument("-ids", help="file with JGI OIDs, one per line. Alternatively can be supplied as a single id or a list of ids separated by '|' at the command prompt", type=str)
    parser.add_argument("-id_type", help="the type of id supplied", choices=['proj', 'img_oid'], required=True)
    parser.add_argument("--get", "-g", help="""folder and regex of data to download separated by '|'. More than one can be included as positional arguments separated by a space Ex: -g 'folder|regex' 'folder|regex'
            \n...Sample: -g 'IMG Data|.*(?:fna|faa|gbk)$' -- matches any file that ends in fna, faa, or gbk in the IMG Data folder.""" , type=arg_type_get, nargs='*')
    parser.add_argument("-out", help="prefix for the downloads", type=str, default=prefix)
    parser.add_argument("-d", "--download", help="single file (as URL) to download")
    parser.add_argument("--r", "--resume", help="resumes file downloads, skips files that already have a local path", action="store_true")
    args = parser.parse_args()

    interface = JGIInterface(args.login, args.ids, args.id_type)

    if args.download:
        interface.download(args.download)
        sys.exit()

    if args.get:
        for genome in interface.genomes:
            genome.download_data(args.get, args.r)
    else:
        interface.genomes[0].print_available_files()

