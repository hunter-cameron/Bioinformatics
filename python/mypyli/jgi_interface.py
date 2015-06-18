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
import datetime
from lxml import etree
import logging

# TODO Add some sort of file check to ensure the entire file was downloaded, specifically for the resume option. 

# TODO Genbank file with > 100 scaffolds -- not sure if this is possible
        # Problem = server-side script checks for scaffold count and demands an email


logging.basicConfig(level=logging.DEBUG)
LOG = logging.getLogger(__name__)

# turn up the level of the requests logger becuase it is noisy
logging.getLogger("requests").setLevel(logging.WARNING)


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
IMG_DATA = 'https://img.jgi.doe.gov/cgi-bin/er/xml.cgi'


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

class AttributeLookupError(AttributeError):
    """ Exception to be thrown with a specific message when a JGIOrganism attribute lookup fails"""
    pass


class JGIBasePath(object):
    """ Base class for files and folders on JGI """

    def __init__(self, name, parent=None):
        self.name = name
        self.parent = parent

        # set the regex match to a full match of the name, if matched later this will be changed
        reg_match = re.match("(.*)", name)



    def get_full_path(self, path_so_far=""):
        """ gets the full path of a piece of a JGI data piece """
        path_so_far = self.name + "/" + path_so_far
        if self.parent is None:
            return "/" + path_so_far
        else:
            return self.get_full_path(path_so_far)

class JGIFolder(JGIBasePath):
    """ Directory on the JGI website, acts as a container for JGI files """


    def __init__(self, name, parent=None):
        super(JGIFolder, self).__init__(name, parent) 
        self.contents = []

    def __str__(self):
        return self.name + "/"

    def append(self, jgi_base_path):
        if isinstance(jgi_base_path, JGIBasePath):
            self.contents.append(jgi_base_path)
        else:
            raise ValueError("jgi_base_path object isn't derrived from JGIBasePath. Refusing to add.")

    def print_tree(self, filehandle=sys.stdout, level=0, recursive=True):
        print("    " * level + str(self))
        level += 1
        for path in self.contents:

            # handle folders
            if isinstance(path, JGIFolder):
                path.print_tree(filehandle, level=level, recursive=True)
            # handle files
            else:
                print("    " * level + str(path))

    def search(self, search_path):
        """ 
        Does a search of a folder and returns all files that match a specified regex.
        Search path can be recursive by supplying a path separated by / such as "myfolder/myfile"
        """
        paths = []
        
        # split off the first portion of the path
        elems = search_path.split("/")


        for path in self.contents:
            match = re.match(elems[0], path.name)

            if match:
                # checkm if this was the last element of the path

                if len(elems) == 1:
                    path.reg_match = match
                    paths.append(path)
                else:
                    if isinstance(path, JGIFolder):
                        try:
                            paths += path.search("/".join(elems[1:]))
                        except DataNotAvailable:
                            pass
        if paths:
            return paths
        else:
            raise DataNotAvailable("No paths found for search '{}'".format(search_path))


class JGIFile(JGIBasePath):
    """ Class to store data about a file from the XML lookup """

    def __init__(self, parent, attribs):
        """ Takes a dict of attributes (from the XML parser) """
        super(JGIFile, self).__init__(attribs["filename"], parent)
        
        self.url = attribs["url"]
        
        # split out the time zone becuase strptime didn't work with it
        ts_elems = attribs["timestamp"].split(" ")
        ts_elems.pop(4)
        ts = " ".join(ts_elems)
        self.timestamp = datetime.datetime.strptime(ts, "%a %b %d %H:%M:%S %Y")
        
        self.size = attribs["size"]
        self.byte_size = int(attribs["sizeInBytes"])
        
        # some entries don't have a md5
        try:
            self.md5 = attribs["md5"]
        except KeyError:
            self.md5 = None

    def __str__(self):
        return "\t".join([self.name, self.size, str(self.timestamp)]) 


class JGIOrganism(object):
    """ Represents an organism entry across JGI and IMG """


    IMG_DATA = ['gbk', 'ko', 'cog', 'pfam', 'tigrfam', 'interpro']


    def __init__(self, interface, **kwargs):
        
        self.interface = interface

        # set download attributes
        prefix = None

        # set JGI attributes
        self.proj_id = None
        self.organism_name = None
        self.jgi_data_tree = None

        # set IMG attributes
        self.taxon_oid = None

         
        valid_kwargs = [
            "prefix", 
            "proj_id", "organism_name",
            "taxon_oid"
            ]
        for k, v in kwargs.items():
            if k not in valid_kwargs:
                raise TypeError("Invalid kwarg: '{}'".format(k))
            else:
                print("Setting self.{} to {}".format(k, v))
                setattr(self, k, v)

    def print_available_files(self, portal="both"):
    
        if portal.lower() == "both":
            # list both
            self._print_img_files()
            print(end="\n\n")
            self._print_jgi_files()
        
        elif portal.lower() == "img":
            # list IMG
            self._print_img_files()
        
        elif portal.lower() == "jgi":
            # list JGI
            self._print_jgi_files()

        else:
            raise ValueError("portal argument must be one of 'both', 'img', or 'jgi'; not '{}'".format(str(portal)))

    def download_data(self, datatype):
        """ 
        Downloads data corresponding to the datatype where datatype is either among the
        elements of IMG_DATA or data is a tuple corresponding to ("path_regex", "file_regex")
        """

        if datatype in self.IMG_DATA:
            # download IMG data
            pass
        
        elif isinstance(datatype, tuple):
            # parse the download regex and then download JGI
            files = self._get_jgi_files(datatype)
            
            try:
                prefix = self.prefix
            except:
                prefix = self.proj_id

            self.interface.download_jgi_files(files=files, prefix=prefix)

        else:
            raise ValueError("Download datatype '{}' is invalid.".format(str(datatype)))
    
    
    
    def _add_children(self, parent, node):
        """
        Recursively adds children nodes to the parent node.

        This could be a classmethod.
        """

        for child in node:
            if child.tag == "file":
                parent.append(JGIFile(parent, child.attrib))

            elif child.tag == "folder":
                folder = JGIFolder(child.attrib["name"], parent)
                parent.append(folder)
                self._add_children(folder, node=child)

            else:
                print("Unknown tag: {}".format(child.tag))

    def _print_jgi_files(self):
        try:
            tree = self.jgi_data_tree
        except AttributeLookupError:
            raise ValueError("jgi_data_tree attribute must be set to print JGI files")

        print("JGI files")
        print("^^^^^^^^^")
        tree.print_tree()
       
    def _print_img_files(self):
        try:
            taxon_oid = self.taxon_oid
        except AttributeLookupError:
            raise ValueError("taxon_oid attrivute must be set to see available IMG files")

        print("IMG files")
        print("^^^^^^^^^")
        for f in self.IMG_DATA:
            print(f)
   
    def _get_jgi_files(self, datatype):
        try:
            taxon_oid = self.proj_id
        except AttributeLookupError:
            raise ValueError("proj_id attribute must be set to download JGI data")

        # split apart the tuple
        folder_reg, file_reg = datatype
        
        # find folders that match the regex
        try:
            folders = self.jgi_data_tree.search(folder_reg)
        except:
            pass

        # find files in the folder
        files= []
        for folder in folders:
            if isinstance(folder, JGIFolder):
                files += folder.search(file_reg)

        return files


    ###################################
    ## Attribute lookups
    #

    def _lookup_proj_id(self):
        """ 
        Returns: the JGI project id as a string.
        Excepts: ValueError
        """
            
        try:
            my_params = {'section': 'TaxonDetail', 'page': 'taxonDetail', 'taxon_oid': self.taxon_oid}
        except AttributeLookupError:
            raise AttributeLookupError("proj_id -> attribute taxon_oid is required to lookup proj_id")

        img_html = self.interface.session.get(IMG_LOOKUP, params=my_params, headers=self.interface.header)
        
        # match the keyValue param in the link with digits to allow the keyValue param to
        # be anywhere within the link 
        # this will need to be changed if JGI ever adds letters to their ids
        m = re.search("\"genome-btn download-btn.*href=['\"].*keyValue=(\d*).*['\"]", img_html.text)
        if m:
            LOG.info("Proj id lookup successful. proj_id={}".format(m.group(1)))
            return m.group(1)
        else:
            #print("Couldn't find URL for IMG id: {}".format(self.taxon_oid), file=sys.stderr)
            raise AttributeLookupError("proj_id -> download link not found for taxon_oid '{}' at URL:\n{}".format(self.taxon_oid, img_html.url))

    def _lookup_organism_name(self):
        """
        Looks up organism name based on project id

        At the name lookup step it will be obvious if there are portal errors.

        """
        
        try:
            my_params = {'keyName': 'jgiProjectId', 'keyValue': self.proj_id}
        except AttributeLookupError:
            raise AttributeLookupError("organism_name -> attribute proj_id must be set to lookup organism name")
        name_html = self.interface.session.get(LOOKUP, params=my_params)    # use the parent session to maintain cookies
        #print((self.id, name_html.url))
       

        # look up the taxon_oid if id isn't already known
        # need some way to handle this if errors happen
        # don't want to 
        if not self.taxon_oid:
            # there are two potential places to get the taxon_oid
            taxon_oid_match = re.search('href="https://img.jgi.doe.gov/genome.php\?id=(\d+)"', name_html.text)
            if taxon_oid_match:
                self.taxon_oid = taxon_oid_match.group(1)
            else:
                taxon_oid_match = re.search('taxon_oid=(\d+)"', name_html.text)

                if taxon_oid_match:
                    self.taxon_oid = taxon_oid_match.group(1)


        # parse the name from the organism info page
        # search for href="/(NAME)/ANYTHING.info.html
        match = re.search('href="/(.*)/.*\.info\.html', name_html.text)

        if match:
            self.name = match.group(1)       # get the first matched group
            return match.group(1)
        else:
            raise AttributeLookupError("organism_name ->  lookup failed for proj_id {}. URL={}".format(self.proj_id, name_html.url))

    def _lookup_taxon_oid(self):
        """
        Looks up taxon_oid based on project id

        """
        try:
            my_params = {'keyName': 'jgiProjectId', 'keyValue': self.proj_id}
        except AttributeLookupError:
            raise AttributeError("taxon_oid -> attribute proj_id must be set to lookup taxon_oid")

        info_html = self.parent.session.get(LOOKUP, params=my_params)    # use the parent session to maintain cookies
       
        # there are two potential places to get the taxon_oid
        taxon_oid_match = re.search('href="https://img.jgi.doe.gov/genome.php\?id=(\d+)"', name_html.text)
        if taxon_oid_match:
            self.taxon_oid = taxon_oid_match.group(1)
        else:
            taxon_oid_match = re.search('taxon_oid=(\d+)"', name_html.text)

            if taxon_oid_match:
                self.taxon_oid = taxon_oid_match.group(1)
            else:
                raise AttributeLookupError("taxon_oid -> lookup failed for project id: {}".format(self.proj_id))

    def _lookup_jgi_data_tree(self):
        """
        Parses the XML structure to return a parent (root) node of a tree of JGIBasePath derrived objects.

        XML Structure Expected:

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
        
        try:
            link = XML + self.organism_name
        except AttributeLookupError:
            raise AttributeLookupError("jgi_data_tree -> attribute organism_name is required for lookup of available files")

        xml = self.interface.session.get(link).text

        root = etree.fromstring(xml)

        # little sanity check to make sure reading the XML format expected
        assert root.tag == "organismDownloads"

        parent = JGIFolder("", parent=None) 

        self._add_children(parent, root)

        return parent 


    @property
    def prefix(self):
        if self.__prefix is None:
            return ""
        else:
            return self.__prefix
    @prefix.setter
    def prefix(self, value):
        self.__prefix = value

    @property
    def proj_id(self):
        if isinstance(self.__proj_id, AttributeLookupError):
            raise self.__proj_id
        
        elif self.__proj_id is None:
            try:
                self.proj_id = self._lookup_proj_id()
            except AttributeLookupError as e:
                self.proj_id = e
                raise
            except RuntimeError:        # this waits for maximum recurson depth, probably a better way
                raise ValueError("Both taxon_oid and proj_id cannot be None. One must be manually set to an appropriate value.")


        
            
        return self.__proj_id
    @proj_id.setter
    def proj_id(self, value):
        self.__proj_id = value
   
    @property
    def organism_name(self):
       
        if isinstance(self.__organism_name, AttributeLookupError):
            raise self.__organism_name

        elif self.__organism_name is None:
            try:
                self.organism_name = self._lookup_organism_name()

            except AttributeLookupError as e:
                self.organism_name = e
                raise
            
        return self.__organism_name
    @organism_name.setter
    def organism_name(self, value):
        self.__organism_name = value

    @property
    def taxon_oid(self):
        if isinstance(self.__taxon_oid, AttributeLookupError):
            raise self.__organism_name

        elif self.__taxon_oid is None:
            print("taxon_oid is NONE")
            try:
                self.taxon_oid = self._lookup_taxon_oid()
            except AttributeLookupError as e:
                self.taxon_oid = e
                raise
            except RuntimeError:        # this waits for maximum recurson depth, probably a better way
                raise ValueError("Both taxon_oid and proj_id cannot be None. One must be manually set to an appropriate value.")

        return self.__taxon_oid
    @taxon_oid.setter
    def taxon_oid(self, value):
        self.__taxon_oid = value

    @property
    def jgi_data_tree(self):
        if isinstance(self.__jgi_data_tree, AttributeLookupError):
            raise self.__jgi_data_tree

        elif self.__jgi_data_tree is None:
            try:
                self.jgi_data_tree = self._lookup_jgi_data_tree()
            except AttributeLookupError as e:
                self.jgi_data_tree = e
                raise

        return self.__jgi_data_tree
    @jgi_data_tree.setter
    def jgi_data_tree(self, value):
        self.__jgi_data_tree = value



class JGIEntry(object):
    """ Class for looking up data on JGI """
    def __init__(self, parent, id, prefix=""):
        self.id = id
        self.parent = parent
        self.prefix = prefix

        # placeholders that will be set during _lookup_name
        self.name = None
        self.data_tree = None
        self.taxon_oid = None




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

        # there are two potential places to get the taxon_oid
        #href="https://img.jgi.doe.gov/genome.php?id=
        taxon_oid_match = re.search('href="https://img.jgi.doe.gov/genome.php\?id=(\d+)"', name_html.text)
        if taxon_oid_match:
            self.taxon_oid = taxon_oid_match.group(1)
        else:
            taxon_oid_match = re.search('taxon_oid=(\d+)"', name_html.text)

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
    def _prepare_entry(self):

        if self.name is None:
            self.name = self._lookup_name()
            if not self.name:
                self.data_tree = None
                return

        self.data_tree = self._parse_xml()

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

        The returned list will be sorted in order of timestamp so it is possible (though a hack)
        to get only most recent file that matches the regex using the -r flag
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
                match = re.match(folder_re, file_name)
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
                new_files.append((directory[file_name], match))
        
        if new_files:
            def _make_datetime(tup):
                
                # split out the time zone
                ts_elems = tup[0]["timestamp"].split(" ")
                ts_elems.pop(4)
                ts = " ".join(ts_elems)

                dt = datetime.datetime.strptime(ts, "%a %b %d %H:%M:%S %Y")
                return dt
            
            # sort new_files by timestamp
            return sorted(new_files, key=_make_datetime, reverse=True)

        else:
            #print("No files found for regex: '{}', id: '{}.'".format(file_reg, self.id), file=sys.stderr)
            #print("Skipping downloading...", file=sys.stderr, end="\n\n")
            #return []

            raise DataNotAvailable("No files found for regex: '{}', id: '{}.'".format(file_reg, self.id))

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

        self._prepare_entry()

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

        if self.taxon_oid:
            return self.taxon_oid
        else:
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
        
        self._prepare_entry() 

        # if there is no tree skip 
        if not self.data_tree:
            print("No files (or access denied) for {}. Skipping Download.".format(self.id))
            return
        
        # add the files to download
        path, file_reg = download_regex
                
        # optionally inject taxon/proj id
        file_reg = file_reg.replace(":taxon_oid", self.lookup_taxon_oid())
        file_reg = file_reg.replace(":proj_id", self.id)
        
        to_download = self._get_download_tuple(path, file_reg)

        # make local path prefix
        if self.prefix:
            prefix = self.prefix
        elif self.name:
            prefix = self.name
        else:
            prefix = self.id

        # make the suffix and download the file
        for attrib, match in to_download:
            url = attrib["url"]
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


    DATA_TYPES = ['fasta', 'gbk', 'ko', 'cog', 'pfam', 'tigrfam', 'interpro']

    def __init__(self, parent, id, prefix=""):
        self.parent = parent
        self.id = id
        self.prefix = prefix
        
        # header to display to the website
        self.header = {'User-Agent': 'Mozilla/5.0 (X11; Linux x86_64; rv:24.0) Gecko/20140924 Firefox/24.0 Iceweasel/24.8.1'}

        #self.name = self._lookup_name()

    def _lookup_proj_id(self):
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
        
        # all of these downloads have a final delim at the end of each line.
        # I should remove that to have properly formatted files
       
        if download_type in self.DATA_TYPES:
            if download_type == 'gbk':
                gbk_file = self.download_genbank()
            elif download_type == 'ko':
                ko_file = self.download_ko()
            elif download_type == 'cog':
                cog_file = self.download_cog()
            elif download_type == 'pfam':
                pfam_file = self.download_pfam()
            elif download_type == 'tigrfam':
                tigrfam_file = self.download_tigrfam()
            elif download_type == 'interpro':
                interpro_file = self.download_interpro()
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
            local_path = self.id + ".gbk"

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
            #print(("pid", pid))
        except:
            #print(request.text)
            
            # attempt to give detailed error mssages
            m = re.search("Please enter your email address since you have selected over 100 entries.", request.text)
            if m:
                raise PortalError("Portal requires email address to download queries with > 100 scaffolds.")

            m = re.search("Please select at least one scaffold", request.text)
            if m:
                #print(request.text)
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

    def download_ko(self):
        """ Method to download a KO file."""
       
       # get the local path
        if self.prefix:
            local_path = self.prefix + ".ko.tsv"
        else:
            local_path = self.id + ".ko.tsv"

        # check the local path
        if not self.parent.check_local_path(local_path):
            return

        print("\nDownloading ko file for {}.".format(self.id), file=sys.stderr)

        # load the KO Ontology page
        payload = {     'section': 'TaxonDetail',
                        'page': 'ko',
                        'taxon_oid': self.id,
                    }
        
        # generate the table and tmpfile 
        request = self.parent.session.post(IMG_LOOKUP, data=payload, headers=self.header)
        
        # parse out the tmpfile
        match = re.search("sid=(yui[^&]*)", request.text)

        try:
            tmpfile = match.group(1)
        except:
            raise PortalError("Error finding KO tmpfile at URL\n{}".format(request.url))
         
        # generate a new payload with the data required to download
        payload = {     'section': 'Selection',
                        'page': 'export',
                        'tmpfile': tmpfile,
                        'table': 'KO',
                        'sort': 'KOID|asc',
                        'rows': 'all',
                        'columns': '[{"key":"Select","label":"Select"},{"key":"KOID","label":"KO ID"},{"key":"Name","label":"Name"},{"key":"Definition","label":"Definition"},{"key":"GeneCount","label":"Gene Count"}]',
                        'c': 'KOID',
                        'f': '',    # left blank on JGI
                        't': 'text',
                        }


        # get the file
        request = self.parent.session.post(IMG_DATA, data=payload, headers=self.header)

        # write it to a file
        with open(local_path, 'w') as OUT:
            OUT.write(request.text)

        print("    Downloaded to {}".format(local_path), file=sys.stderr)

    def download_cog(self):
        """ Method to download a COG summary file."""
       
       # get the local path
        if self.prefix:
            local_path = self.prefix + ".cog.tsv"
        else:
            local_path = self.id + ".cog.tsv"

        # check the local path
        if not self.parent.check_local_path(local_path):
            return

        print("\nDownloading cog file for {}.".format(self.id), file=sys.stderr)

        # load the COG page
        payload = {     'section': 'TaxonDetail',
                        'page': 'cogs',
                        'taxon_oid': self.id,
                    }
        
        # generate the table and tmpfile 
        request = self.parent.session.post(IMG_LOOKUP, data=payload, headers=self.header)
        
        # parse out the tmpfile
        match = re.search("sid=(yui[^&]*)", request.text)

        try:
            tmpfile = match.group(1)
        except:
            raise PortalError("Error finding COG tmpfile at URL\n{}".format(request.url))
         
        # generate a new payload with the data required to download
        payload = {     'section': 'Selection',
                        'page': 'export',
                        'tmpfile': tmpfile,
                        'table': 'cog',
                        'sort': 'COGID|asc',
                        'rows': 'all',
                        'columns': '[{"key":"Select","label":"Select"},{"key":"COGID","label":"COG ID"},{"key":"COGName","label":"COG Name"},{"key":"GeneCount","label":"Gene Count"}]',
                        'c': 'COGID',
                        'f': '',    # left blank on JGI
                        't': 'text',
                        }


        # get the file
        request = self.parent.session.post(IMG_DATA, data=payload, headers=self.header)

        # write it to a file
        with open(local_path, 'w') as OUT:
            OUT.write(request.text)

        print("    Downloaded to {}".format(local_path), file=sys.stderr)

    def download_pfam(self):
        """ Method to download a Pfam summary file."""
       
       # get the local path
        if self.prefix:
            local_path = self.prefix + ".pfam.tsv"
        else:
            local_path = self.id + ".pfam.tsv"

        # check the local path
        if not self.parent.check_local_path(local_path):
            return

        print("\nDownloading pfam file for {}.".format(self.id), file=sys.stderr)

        # load the pfam  page
        payload = {     'section': 'TaxonDetail',
                        'page': 'pfam',
                        'taxon_oid': self.id,
                    }
        
        # generate the table and tmpfile 
        request = self.parent.session.post(IMG_LOOKUP, data=payload, headers=self.header)
        
        # parse the sid from the tmpfile
        match = re.search("sid=(yui[^&]*)", request.text)

        try:
            tmpfile = match.group(1)
        except:
            raise PortalError("Error finding Pfam tmpfile at URL\n{}".format(request.url))
         
        # generate a new payload with the data required to download
        payload = {     'section': 'Selection',
                        'page': 'export',
                        'tmpfile': tmpfile,
                        'table': 'pfam',
                        'sort': 'PfamID|asc',
                        'rows': 'all',
                        'columns': '[{"key":"Select","label":"Select"},{"key":"PfamID","label":"Pfam ID"},{"key":"PfamName","label":"Pfam Name"},{"key":"GeneCount","label":"Gene Count"}]',
                        'c': 'PfamID',
                        'f': '',    # left blank on JGI
                        't': 'text',
                        }


        # get the file
        request = self.parent.session.post(IMG_DATA, data=payload, headers=self.header)

        # write it to a file
        with open(local_path, 'w') as OUT:
            OUT.write(request.text)

        print("    Downloaded to {}".format(local_path), file=sys.stderr)

    def download_tigrfam(self):
        """ Method to download a TIGRfam summary file."""
       
       # get the local path
        if self.prefix:
            local_path = self.prefix + ".tigrfam.tsv"
        else:
            local_path = self.id + ".tigrfam.tsv"

        # check the local path
        if not self.parent.check_local_path(local_path):
            return

        print("\nDownloading tigrfam file for {}.".format(self.id), file=sys.stderr)

        # load the tigrfam page
        payload = {     'section': 'TaxonDetail',
                        'page': 'tigrfam',
                        'taxon_oid': self.id,
                    }
        
        # generate the table and tmpfile 
        request = self.parent.session.post(IMG_LOOKUP, data=payload, headers=self.header)
        
        # parse the sid from the tmpfile
        match = re.search("sid=(yui[^&]*)", request.text)

        try:
            tmpfile = match.group(1)
        except:
            raise PortalError("Error finding tigrfam tmpfile at URL\n{}".format(request.url))
         
        # generate a new payload with the data required to download
        payload = {     'section': 'Selection',
                        'page': 'export',
                        'tmpfile': tmpfile,
                        'table': 'tigrfam',
                        'sort': 'TIGRfamID|asc',
                        'rows': 'all',
                        'columns': '[{"key":"Select","label":"Select"},{"key":"TIGRfamID","label":"TIGRfam ID"},{"key":"TIGRfamName","label":"TIGRfam Name"},{"key":"GeneCount","label":"Gene Count"}]',
                        'c': 'TIGRfamID',
                        'f': '',    # left blank on JGI
                        't': 'text',
                        }


        # get the file
        request = self.parent.session.post(IMG_DATA, data=payload, headers=self.header)

        # write it to a file
        with open(local_path, 'w') as OUT:
            OUT.write(request.text)

        print("    Downloaded to {}".format(local_path), file=sys.stderr)

    def download_interpro(self):
        """ Method to download an InterPro summary file."""
       
       # get the local path
        if self.prefix:
            local_path = self.prefix + ".interpro.tsv"
        else:
            local_path = self.id + ".interpro.tsv"

        # check the local path
        if not self.parent.check_local_path(local_path):
            return

        print("\nDownloading interpro file for {}.".format(self.id), file=sys.stderr)

        # load the InterPro page
        payload = {     'section': 'TaxonDetail',
                        'page': 'ipr',
                        'taxon_oid': self.id,
                    }
        
        # generate the table and tmpfile 
        request = self.parent.session.post(IMG_LOOKUP, data=payload, headers=self.header)
        
        # parse the sid from the tmpfile
        match = re.search("sid=(yui[^&]*)", request.text)

        try:
            tmpfile = match.group(1)
        except:
            raise PortalError("Error finding InterPro tmpfile at URL\n{}".format(request.url))
         
        # generate a new payload with the data required to download
        payload = {     'section': 'Selection',
                        'page': 'export',
                        'tmpfile': tmpfile,
                        'table': 'interpro',
                        'sort': 'InterProID|asc',
                        'rows': 'all',
                        'columns': '[{"key":"Select","label":"Select"},{"key":"InterProID","label":"InterPro ID"},{"key":"InterProName","label":"InterPro Name"},{"key":"GeneCount","label":"Gene Count"}]',
                        'c': 'TIGRfamID',
                        'f': '',    # left blank on JGI
                        't': 'text',
                        }


        # get the file
        request = self.parent.session.post(IMG_DATA, data=payload, headers=self.header)

        # write it to a file
        with open(local_path, 'w') as OUT:
            OUT.write(request.text)

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

        # header to display to the website
        self.header = {'User-Agent': 'Mozilla/5.0 (X11; Linux x86_64; rv:24.0) Gecko/20140924 Firefox/24.0 Iceweasel/24.8.1'}




    @staticmethod
    def _login(login_file=''):
        """
        Parses user credentials from a file or string and then uses them to set up a session. Returns the session,
        otherwise dies.

        """

        #parse the login info
        if login_file == '':
            raise ValueError("Login file not specified")


        # allow this method to take a username/pass as a string
        try:
            with open(login_file, 'r') as in_handle:
                login_data = {
                    'login': in_handle.readline().strip(),
                    'password':  in_handle.readline().strip()
                }

        except FileError:
            login_elems = login_file.split("\n")
            login_data = {
                    'login': login_elems[0].strip(),
                    'password': login_elems[1].strip()
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
        except (FileNotFoundError, OSError):
            arg = arg.strip()
            items = arg.split(delim)
            
        return items

    @staticmethod
    def _resolve_path_conflicts(files, newest_only=False):
        """ 
        Resolves multiple files attempting to be downloaded to the same local path.
        Either adds the file date or keeps only the newest
        """

        # map file to suffix
        file_dict = {}
        for file in files:
            try:
                suffix = file.reg_match.group(1)
            except IndexError:
                suffix = file.name
            file_dict[file] = suffix

        # find any repeated suffices
        suf_count = {}
        for suf in file_dict.values():
            suf_count[suf] = suf_count.get(suf, 0) + 1
        repeated_sufs = [suf for suf in suf_count if suf_count[suf] > 1]


        #
        ## Make a new dict with unique sufs
        #

        # add the unique sufs
        # this line is arbitrarily complex. Basically just transfers any key with a value that is
        # not in the repeated_sufs list to the new dict
        new_f_dict = {k: file_dict[k] for k in file_dict if file_dict[k] not in repeated_sufs}
        
        # handle any duplicated sufs 
        for suf in [suf for suf in suf_count if suf_count[suf] > 1]:
            files_w_suf = [f for f in file_dict if file_dict[f] == suf]
            
            if newest_only:
                newest = sorted(files_w_suf, key=lambda x: x.timestamp, reverse=True)[0]
                new_f_dict[newest] = suf
            else:
                # prepend date to the suffix
                for file in files_w_suf:
                    new_f_dict[file] = file.timestamp.strftime("%Y%b%d") + "." + suf
                
            
        #### Check that new_f_dict is unique now

        return new_f_dict

    def download_jgi_files(cls, files, prefix="organism"):
        file_dict = cls._resolve_path_conflicts(files, newest_only=True)
        
        for file, suffix in file_dict.items():
            local_path = prefix + "." + suffix
            if cls.check_local_path(local_path):
                cls._download_data(DATA + file.url, local_path)
        

    def _download_data(self, url, local_path):
        """ 
        Downloads a file in chunks.
        Returns: Nothing
        Excepts: AssertionError
        """

        #download the file in chunks

        if self.check_local_path(local_path):

            print("Downloading {} to\n    {}".format(url, local_path), end="\n\n")
            request = self.session.get(url, stream=True)
            
            try:
                with open(local_path, 'wb') as fh:
                    for chunk in request.iter_content(chunk_size=1024):
                        if chunk:  # filter out keep-alive new chunks
                            fh.write(chunk)
                            fh.flush()
            except IOError:
                print("Problem downloading to local path {}".format(local_path), file=sys.stderr)

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
                        try:
                            download_entry = interface.convert_entry(entry)
                        except PortalError as e:
                            print(e, file=sys.stderr)
                            failed_log.append(download_regex)
                            continue
                    else:
                        raise ValueError("Incorrect --get type for your id. Specify --convert to convert between ids.")

                # try to download data
                try:
                    download_entry.download_data(download_regex)
            
                except (DataNotAvailable, PortalError) as e:
                    print(e, file=sys.stderr, end="\n\n")
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
        failed_ids = []
        for entry in interface.entries:
            print("Printing available files for all ids. Press Ctrl+c to interrupt.", file=sys.stderr)
            entry.print_available_files()

            if args.convert:
                try:
                    conv_entry = interface.convert_entry(entry)
                except PortalError as e:
                    print(e, file=sys.stderr, end="\n\n")
                    failed_ids.append(entry)
                    continue
                conv_entry.print_available_files()

        if failed_ids:
            print("\n\nThere were errors in fetching the available files for the following ids:")
            [print("    " + entry.id) for entry in failed_ids]


def main2(args):
    """
    Main function to manipulate the classes to download data.
    
    Maybe should be incorporated into JGIInterface.
    """
    
    interface = JGIInterface(args.login, force=args.force, convert=args.convert, resume=args.resume)
    
    org = JGIOrganism(interface, taxon_oid=args.ids)
    return org
    print(org.taxon_oid)
    print(org.proj_id)   
    sys.exit()


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

    org = main2(args)

