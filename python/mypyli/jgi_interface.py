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
from html.parser import HTMLParser
import traceback

# TODO Add some sort of file check to ensure the entire file was downloaded, specifically for the resume option. 

# TODO Genbank file with > 100 scaffolds -- not sure if this is possible
        # Problem = server-side script checks for scaffold count and demands an email

# TODO With JGI's initiative to make portals for all IMG data, there are some weird naming conventions for things that were sequenced at other locations so do not have a JGI portal. So best I can tell, such projects do not have a project name but the download link (with a different format) gives the organism lookup name directly.
        # Right now, this isn't a big enough problem to address but, if it becomes bigger, something will need to be done.

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
IMG_LOOKUP = 'https://img.jgi.doe.gov/cgi-bin/mer/main.cgi'
#IMG_LOOKUP = 'https://img.jgi.doe.gov/cgi-bin/er/main.cgi'
IMG_DATA = 'https://img.jgi.doe.gov/cgi-bin/mer/xml.cgi'


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
            return path_so_far
        else:
            return self.parent.get_full_path(path_so_far)

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
    """ Stores data about a file from the XML lookup """

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

class MetadataParser(HTMLParser):
    def __init__(self):
        super(MetadataParser, self).__init__()
        self.data_dict = {}

        # stores the current data key
        self.cur_key = None

        # stores what to do with the data
        self.slurp_mode = None

    def handle_starttag(self, tag, attr):
        if tag == "th":
            # headers are the keys
            self.slurp_mode = "key"
        elif tag == "td":
            self.slurp_mode = "data"

    def handle_data(self, data):
        data = data.strip()
        if data:
            if self.slurp_mode == "key":
                self.cur_key = data.strip()
            elif self.slurp_mode == "data":
                self.data_dict[self.cur_key] = self.data_dict.get(self.cur_key, "") + data.strip()

    def get_metadata(self, table_html):
        self.feed(table_html)
        
        # map JGI's names to more friendly ones
        rough_to_standard = {
                "Study Name (Proposal Name)": "proposal",
                "Organism Name": "organism_name",
                "Taxon ID": "taxon_oid",
                "IMG Submission ID": "submission_id",
                "NCBI Taxon ID": "ncbi_tax_id",
        #        "GOLD ID in IMG Database": "",
        #        "GOLD Analysis Project Id": "",
        #        "GOLD Analysis Project Type": "",
                "Submission Type": "submission",
        #        "External Links": "",
                "Lineage": "lineage",
                "Sequencing Status": "seq_status",
                "Sequencing Center": "seq_center",
        #        "IMG Release": "",
        #        "Comment": "",
        #        "Release Date": "",
        #        "Add Date": "",
        #        "Modified Date": "",
        #        "Distance Matrix Calc. Date": "",
                "High Quality": "high_quality",
        #        "IMG Product Flag": "img_product",
                "Is Public": "public",
        #        "Project Information": "",
        #        "Bioproject Accession": "",
        #        "Biosample Accession": "",
                "Culture Type": "culture_type",
        #        "Uncultured Type": "culture_type",     # for single cells, their type is listed as uncultured
                "GOLD Sequencing Strategy": "seq_strategy",
                "Gram Staining": "gram",
        #        "Seq Status": "",
                "Sequencing Method": "seq_method",
        #        "Type Strain": ""
            }


        # instead of using nice names, maybe it would be best to use the full JGI name?
        # that would be most conducive to essentially having a local copy of JGI's data
        # I still only want a subset of the headers
        interesting_data = {
                "Study Name (Proposal Name)",
                "Organism Name",
                "Taxon ID",
                "IMG Submission ID",
                "NCBI Taxon ID",
        #        "GOLD ID in IMG Database",
        #        "GOLD Analysis Project Id",
        #        "GOLD Analysis Project Type",
                "Submission Type",
        #        "External Links",
                "Lineage",
                "Sequencing Status",
                "Sequencing Center",
        #        "IMG Release",
        #        "Comment",
        #        "Release Date",
        #        "Add Date",
        #        "Modified Date",
        #        "Distance Matrix Calc. Date",
                "High Quality",
        #        "IMG Product Flag",
                "Is Public",
        #        "Project Information",
        #        "Bioproject Accession",
        #        "Biosample Accession",
                "Culture Type",
                "Uncultured Type",      # single cells have this field instead
                "GOLD Sequencing Strategy",
                "Gram Staining",
        #        "Seq Status",
                "Sequencing Method",
        #        "Type Strain",
            }


        # make a polished dict (to give some consistency to available data and reduce KeyErrors)
        polished_dict = {}
        for key in interesting_data:
            polished_dict[key] = self.data_dict.get(key, None)

        return polished_dict

class JGIOrganism(object):
    """ 
    Represents an organism entry across JGI and IMG
    
    This class handles the entry until the data has been chosen for download.
    Then, the data is passed to the JGIInterface for local download. 
    """


    IMG_DATA = ['gbk', 'ko', 'cog', 'pfam', 'tigrfam', 'interpro']
    IMG_DATA_SUF = {
            "gbk": "gbk",
            "ko": "ko.txt",
            "cog": "cog.txt", 
            "pfam": "pfam.txt",
            "tigrfam": "tigrfam.txt",
            "interpro": "interpro.txt"
            }


    def __init__(self, interface, **kwargs):
        
        self.interface = interface

        # reports to store to give information about access
        self.download_reports = {}

        # set download attributes
        self.prefix = None

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
                LOG.debug("Setting self.{} to {}".format(k, v))
                setattr(self, k, v)

    def __str__(self):
        """
        This method must tread carefully, I don't want to trigger and data lookups.
        I only want to print what data is available.
        """
        string = "JGIOrganism Object:"

        # if the user has given a prefix to this Organism, just use that
        if self.prefix:
            return string + " " + self.prefix

        for attr in ["taxon_oid", "proj_id", "organism_name"]:
            try:
                value = getattr(self, "_" + attr)
                # don't want to print it if it is an error
            except AttributeError:
                continue
            if isinstance(value, str):
                string += "\t{}:{}".format(attr, value)
            else:
                continue
        return string
   
    def print_download_report(self):
        for k in sorted(self.download_reports):
            print(k, end="")
            
            if isinstance(self.download_reports[k], dict):
                print()
                for k2 in sorted(self.download_reports[k]):
                    print("\t{}\t{}".format(k2, self.download_reports[k][k2]))
            else:
                print("\t{}".format(self.download_reports[k]))
            
    def print_available_files(self, portal="both"):
    
        print("####################", end="\n\n")
        if portal.lower() == "both":
            # list both
            self._print_img_files()
            self._print_jgi_files()
        
        elif portal.lower() == "img":
            # list IMG
            self._print_img_files()
        
        elif portal.lower() == "jgi":
            # list JGI
            self._print_jgi_files()

        else:
            raise ValueError("portal argument must be one of 'both', 'img', or 'jgi'; not '{}'".format(str(portal)))
        print("####################")

    def get_metadata(self):
        """ Gets a bunch of metadata from the IMG portal page. Returns is as a dict """

        try:
            my_params = {'section': 'TaxonDetail', 'page': 'taxonDetail', 'taxon_oid': self.taxon_oid}
        except AttributeLookupError:
            raise AttributeLookupError("proj_id -> attribute taxon_oid is required to lookup proj_id")

        LOG.debug("Getting metadata for {}".format(str(self)))

        img_portal = self.interface.session.get(IMG_LOOKUP, params=my_params, headers=self.interface.header)
        
        # get just the data table
        match = re.search("<p></p><a name='overview' href='#'><h2>Overview</h2> </a>(.*?)</table>", img_portal.text, re.DOTALL)
        if match:
            md_table = match.group(1)
        else:
            raise PortalError("Could not find the metadata table. The HTML is different than expected (possibly because the lookup failed.\n{}".format(str(self)))

        parser = MetadataParser()
        return parser.get_metadata(md_table)


    def download_data(self, datatype):
        """ 
        Downloads data corresponding to the datatype where datatype is either among the
        elements of IMG_DATA or data is a tuple corresponding to ("path_regex", "file_regex")
        """
 
        LOG.debug("Trying to download '{}' for {}".format(str(datatype), str(self))) 
        if datatype in self.IMG_DATA:
            try:
                # download IMG data
                tmpfile = self._img_data_request(datatype)
                payload = self._get_img_download_payload(datatype, tmpfile)
            except (AttributeLookupError, PortalError) as e:
                self.download_reports[datatype] = "failed"
                raise 

            # set the prefix
            try:
                prefix = self.prefix
            except:
                prefix = self.taxon_oid

            self.download_reports[datatype] = self.interface.download_img_file(payload, prefix, datatype)
            

        elif isinstance(datatype, tuple):
            try:
            # parse the download regex and then download JGI
                files = self._get_jgi_files(datatype)
            except (AttributeLookupError, PortalError) as e:
                self.download_reports["|".join(datatype)] = "failed"
            try:
                prefix = self.prefix
            except:
                prefix = self.proj_id


            self.download_reports["|".join(datatype)] = self.interface.download_jgi_files(files=files, prefix=prefix)

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

        LOG.info("Printing available JGI files for '{}'".format(str(self)))
        print("JGI files")
        print("^^^^^^^^^")
        tree.print_tree()
       
    def _print_img_files(self):
        try:
            taxon_oid = self.taxon_oid
        except AttributeLookupError:
            raise ValueError("taxon_oid attrivute must be set to see available IMG files")

        LOG.info("Printing available IMG files for '{}'".format(str(self)))
        print("IMG files")
        print("^^^^^^^^^")
        for f in self.IMG_DATA:
            print(f)
   
    def _get_jgi_files(self, datatype):
        LOG.debug("Getting files from JGI that match the regex '{}'".format(str(datatype)))
        try:
            pass
            #taxon_oid = self.proj_id
        except AttributeLookupError:
            raise ValueError("proj_id attribute must be set to download JGI data")

        # split apart the tuple
        folder_reg, file_reg = datatype
        
        # find folders that match the regex
        try:
            folders = self.jgi_data_tree.search(folder_reg)
        except DataNotAvailable:
            raise DataNotAvailable("No folder matching '{}' was found.".format(folder_reg))

        # find files in the folder
        files= []
        for folder in folders:
            if isinstance(folder, JGIFolder):
                try:
                    files += folder.search(file_reg)
                except DataNotAvailable:
                    raise DataNotAvailable("File matching '{}' was not found in folder '{}'".format(file_reg, folder.get_full_path()))

        

        return files

    def _get_img_files(self, datatype):

        try:
            taxon_oid = self.taxon_oid
        except AttributeLookupError:
            raise ValueError("attribute taxon_oid must be set to download IMG data")

        if datatype == "gbk":
            pass
        elif datatype in ["ko", "cog", "pfam", "tigrfam", "interpro"]:
            payloads = self._get_img_payloads(datatype)

    def _img_data_request(self, datatype):
        """ Makes an IMG data request and returns a pid that can be used to access the data on the server """
        taxon_oid = self.taxon_oid

        # generate a payload to ask the server for a data file
        if datatype == "gbk":
            payload = {     'taxon_oid': taxon_oid,
                            'scaffold_oid': 'all',
                            'format': 'gbk',
                            '_section_TaxonDetail_processArtemisFile': 'Go'
                        }

            pid_reg = "id='pid' name='pid' value=\'(\d+)\'"

        elif datatype in ['ko', 'cog', 'pfam', 'tigrfam', 'interpro']:
            payload = {     'section': 'TaxonDetail',
                            'taxon_oid': taxon_oid,
                        }
        
            page_dict = {
                    'ko': 'ko',
                    'cog': 'cogs',
                    'pfam': 'pfam',
                    'tigrfam': 'tigrfam',
                    'interpro': 'ipr'
                    }

            payload['page'] = page_dict[datatype]
        
            pid_reg = "sid=(yui[^&]*)"
       

        # request the data be generated
        request = self.interface.session.post(IMG_LOOKUP, data=payload, headers=self.interface.header)
 
        # get the pid
        match = re.search(pid_reg, request.text)

        try:
            pid = match.group(1)
            return pid
        except:
            #print(request.text)
            # attempt to give detailed error mssages

            # some genbank specific messages
            m = re.search("Please enter your email address since you have selected over 100 entries.", request.text)
            if m:
                raise PortalError("Portal requires email address to download queries with > 100 scaffolds.")

            m = re.search("Please select at least one scaffold", request.text)
            if m:
                #print(request.text)
                raise PortalError("Portal claims no scaffolds were selected (this is a true portal error?).")
            
            
            LOG.debug(request.url)
            print(request.text)
            raise PortalError("PID not found for datatype '{}' in {}.".format(datatype, str(self)))

    @staticmethod
    def _get_img_download_payload(datatype, tmpfile):
        """ Looks up the payloads for various IMG data types. """
       
        if datatype in ["ko", "cog", "pfam", "tigrfam", "interpro"]:
            payload = { 

                        'section': 'Selection',
                        'page': 'export',
                        'tmpfile': tmpfile,
                        'rows': 'all',
                        'f': '',    # left blank on JGI
                        't': 'text',
                    }
        else:
            payload = {}


        additional_data = {

                "gbk": {
                        'pid': tmpfile,
                        '_section_TaxonDetail_downloadArtemisFile_noHeader': 'Download File',
                        'type': 'gbk'
                    },

                        

                "ko": {
                        'table': 'KO',
                        'sort': 'KOID|asc',
                        'columns': '[{"key":"Select","label":"Select"},{"key":"KOID","label":"KO ID"},{"key":"Name","label":"Name"},{"key":"Definition","label":"Definition"},{"key":"GeneCount","label":"Gene Count"}]',
                        'c': 'KOID'
                    },


                "cog": {
                        'table': 'cog',
                        'sort': 'COGID|asc',
                        'rows': 'all',
                        'columns': '[{"key":"Select","label":"Select"},{"key":"COGID","label":"COG ID"},{"key":"COGName","label":"COG Name"},{"key":"GeneCount","label":"Gene Count"}]',
                        'c': 'COGID'
                    },

                    "pfam": {
                        'table': 'pfam',
                        'sort': 'PfamID|asc',
                        'rows': 'all',
                        'columns': '[{"key":"Select","label":"Select"},{"key":"PfamID","label":"Pfam ID"},{"key":"PfamName","label":"Pfam Name"},{"key":"GeneCount","label":"Gene Count"}]',
                        'c': 'PfamID',
                    },

                "tigrfam": {
                        'table': 'tigrfam',
                        'sort': 'TIGRfamID|asc',
                        'rows': 'all',
                        'columns': '[{"key":"Select","label":"Select"},{"key":"TIGRfamID","label":"TIGRfam ID"},{"key":"TIGRfamName","label":"TIGRfam Name"},{"key":"GeneCount","label":"Gene Count"}]',
                        'c': 'TIGRfamID',
                    },

                "interpro": {
                        'table': 'interpro',
                        'sort': 'InterProID|asc',
                        'rows': 'all',
                        'columns': '[{"key":"Select","label":"Select"},{"key":"InterProID","label":"InterPro ID"},{"key":"InterProName","label":"InterPro Name"},{"key":"GeneCount","label":"Gene Count"}]',
                        'c': 'InterProID',
                    },
            }
        

        payload.update(additional_data[datatype])

        return payload


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
            LOG.debug("Proj id lookup successful. proj_id={}".format(m.group(1)))
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

        info_html = self.interface.session.get(LOOKUP, params=my_params)    # use the parent session to maintain cookies
       
        # there are two potential places to get the taxon_oid
        taxon_oid_match = re.search('href="https://img.jgi.doe.gov/genome.php\?id=(\d+)"', info_html.text)
        if taxon_oid_match:
            return taxon_oid_match.group(1)
        else:
            taxon_oid_match = re.search('taxon_oid=(\d+)"', info_html.text)

            if taxon_oid_match:
                return taxon_oid_match.group(1)
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
        if self._prefix is None:
            return ""
        else:
            return self._prefix
    @prefix.setter
    def prefix(self, value):
        self._prefix = value

    @property
    def proj_id(self):
        if isinstance(self._proj_id, AttributeLookupError):
            raise self._proj_id
        
        elif self._proj_id is None:
            try:
                self.proj_id = self._lookup_proj_id()
            except AttributeLookupError as e:
                self.proj_id = e
                raise
            except RuntimeError:        # this waits for maximum recurson depth, probably a better way
                raise ValueError("Both taxon_oid and proj_id cannot be None. One must be manually set to an appropriate value.")

            
        return self._proj_id
    @proj_id.setter
    def proj_id(self, value):
       self._proj_id = value
   
    @property
    def organism_name(self):
       
        if isinstance(self._organism_name, AttributeLookupError):
            raise self._organism_name

        elif self._organism_name is None:
            try:
                self.organism_name = self._lookup_organism_name()

            except AttributeLookupError as e:
                self.organism_name = e
                raise
            
        return self._organism_name
    @organism_name.setter
    def organism_name(self, value):
        self._organism_name = value

    @property
    def taxon_oid(self):
        if isinstance(self._taxon_oid, AttributeLookupError):
            raise self._taxon_oid

        elif self._taxon_oid is None:
            try:
                self.taxon_oid = self._lookup_taxon_oid()

            except AttributeLookupError as e:
                self.taxon_oid = e
                raise
            except RuntimeError:        # this waits for maximum recurson depth, probably a better way
                raise ValueError("Both taxon_oid and proj_id cannot be None. One must be manually set to an appropriate value.")

        return self._taxon_oid
    @taxon_oid.setter
    def taxon_oid(self, value):
        self._taxon_oid = value

    @property
    def jgi_data_tree(self):
        if isinstance(self._jgi_data_tree, AttributeLookupError):
            raise self._jgi_data_tree

        elif self._jgi_data_tree is None:
            try:
                self.jgi_data_tree = self._lookup_jgi_data_tree()
            except AttributeLookupError as e:
                self.jgi_data_tree = e
                raise

        return self._jgi_data_tree
    @jgi_data_tree.setter
    def jgi_data_tree(self, value):
        self._jgi_data_tree = value



class JGIInterface(object):
    """
    Manages the connection with JGI, responsible for logging in, 
    establishing, and maintaining cookies throughout the session
    """

    def __init__(self, **kwargs):
        
        # need a method of reporting download history (successful, failed, etc)
        self.failed_downloads = []
        self.overwritten_paths = []
        self.skipped_paths = []

        # set default values for kwargs
        self.force_overwrite = False
        self.resume = False
        self.newest_only = False

        # header to display to the website
        self.header = {'User-Agent': 'Mozilla/5.0 (X11; Linux x86_64; rv:24.0) Gecko/20140924 Firefox/24.0 Iceweasel/24.8.1'}

        # handle kwargs
        valid_kwargs = [
            "force_overwrite",
            "resume",
            "newest_only",
            "login_file",
            "username",
            "password"
            ]
        for k, v in kwargs.items():
            if k not in valid_kwargs:
                raise TypeError("Invalid kwarg: '{}'".format(k))
            else:
                if k in ["force_overwrite", "resume", "newest_only"]:
                    if isinstance(v, bool):
                        LOG.debug("Setting interface.{} to {}".format(k, v))
                        setattr(self, k, v)
                    else:
                        raise ValueError("kwarg '{}' must be a boolean".format(k))
                if k == "login_file":
                    with open(v, 'r') as IN:
                        username = IN.readline().strip()
                        password = IN.readline().strip()
                    
                    self.session = self._login(username, password)

                if k == "username":
                    try:
                        self.session = self._login(v, kwargs["password"])
                    except KeyError:
                        raise ValueError("when kwarg 'username' is given kwarg 'password' must also be given")

                    
                #setattr(self, k, v)

        # make sure session was set
        try:
            self.session
        except AttributeError:
            raise ValueError("JGI login information is required. Please specify kwargs 'username' and 'password' OR kwarg 'login_file'")

    @staticmethod
    def _login(username, password):
        """
        Parses user credentials from a file or string and then uses them to set up a session. Returns the session,
        otherwise dies.

        """

        LOG.info("Logging into JGI...")
        # convert user and pass to a payload
        login_data = {
            'login': username,
            'password':  password
            }

        # make a session to store cookie information like a browser would
        session = requests.session()

        # login
        request = session.post(LOGIN, data=login_data)

        # check status of logon
        match = re.search('You have signed in successfully.', request.text)
        if not match:
            print(request.text)       # uncomment to display raw login info. Useful to look through if having problems
            raise PortalError("Login may not have been successful.")
        else:
            LOG.info("Login successful!")
            return session

    def get_taxon_oids_for_proposal(self, proposal_name):
        """ Returns a list of JGIOrganisms created from the IMG search results of a proposal name. """


        """
        On the IMG site, you must search in one action and then add extra fields in another. 
        I would like to do this in one action.

        However, the field addition also passes two keys with the same name. Perhaps I can pass an array?
        Yes, I can pass array but then I must post as JSON.
        I can also use a list of tuples for the params instead of dict and use data insead of json.
        """
        
        payload = [ ('section', 'FindGenomes'), 
                    ('page', 'findGenomeResults'),
                    ('taxonSearchTerm', proposal_name),
                    ('taxonSearchFilter', 'proposal_name'),
                    ('find_organism', 'Go'),
                    
                    # These two commands are part of the second wave that I thought I might 
                    # could just add to the first, didn't work. So, for now, I'll omit it
                    #('genome_field_col', 't.taxon_oid'),
                    #('genome_field_col', 't.jgi_project_id'),
                    ]

        request = self.session.post(IMG_LOOKUP, data=payload, headers=self.header)

        # parse temporary id
        pid_reg = "sid=(yui[^&]*)"
        match = re.search(pid_reg, request.text)

        try:
            tmpid = match.group(1)
        except:
            raise JGIPortalError("Temporary id not found while looking up proposals for '{}'".format(proposal_name))

        download_payload = {
                    'section': 'Selection',
                    'page': 'export',
                    'tmpfile': tmpid,
                    'table': 'taxontable',
                    'sort': 'domain|asc',
                    'rows': 'all',

                    # original
                    #'columns': '[{"key":"Select","label":"Select"},{"key":"Domain","label":"Domain"},{"key":"Status","label":"Status"},{"key":"StudyName","label":"Study Name"},{"key":"GenomeNameSampleName","label":"Genome Name / Sample Name"},{"key":"SequencingCenter","label":"Sequencing Center"},{"key":"IMGGenomeID","label":"IMG Genome ID "},{"key":"JGIProjectID","label":"JGI Project ID"},{"key":"GenomeSize","label":"Genome Size "},{"key":"GeneCount","label":"Gene Count "}]',

                    # truncated to only taxon_oid -- this was done because most of the fields
                    # in the column list above were blank (and all I NEED is taxon_oid)
                    'columns': '[{"key":"Select","label":"Select"}]',
                    'c': 'domain',
                    'f': '',
                    't': 'text',
                    }
        
        request = self.session.post(IMG_DATA, data=download_payload, headers=self.header)

        lines = request.text.split("\n")
        
        # make a list of taxon_oids skipping the header line and any blank lines
        taxon_oids = []
        for taxon_oid in lines[1:]:
            taxon_oid = taxon_oid.strip()
            if taxon_oid:
                taxon_oids.append(taxon_oid)

        return taxon_oids


    def _resolve_path_conflicts(self, files):
        """ 
        Resolves multiple files attempting to be downloaded to the same local path.
        Either adds the file date or keeps only the newest.
        """

        LOG.debug("Checking for and resolving path conflicts...")
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
            
            if self.newest_only:
                newest = sorted(files_w_suf, key=lambda x: x.timestamp, reverse=True)[0]
                new_f_dict[newest] = suf
            else:
                # prepend date to the suffix
                for file in files_w_suf:
                    new_f_dict[file] = file.timestamp.strftime("%Y%b%d") + "." + suf
                
            
        #### Check that new_f_dict is unique now

        return new_f_dict

    def download_jgi_files(self, files, prefix="organism"):
        """ 
        Accepts a list of JGIFile objects.
        Checks file naming for local collision among the batch and resolves conflicts.
        Then, forwards each file to _download_data
        """

        # begin making a download report
        download_report = {str(file): "pending" for file in files}

        LOG.debug("Downloading JGI files...")
        file_dict = self._resolve_path_conflicts(files)
        
        # report how many files were removed due to name conflict
        for file in download_report:
            if str(file) not in file_dict:
                download_report[file] = "name conflict"

        for file, suffix in file_dict.items():
            local_path = prefix + "." + suffix
            if self._proceed_with_download(local_path, file=file):
                self._download_data(DATA + file.url, local_path)
                download_report[str(file)] = "successful"
            else:
                download_report[str(file)] = "skipped" 

        return download_report

    def download_img_file(self, payload, prefix, datatype):

        local_path = prefix + "." + JGIOrganism.IMG_DATA_SUF[datatype]
        if datatype == "gbk":
            request = self.session.post(IMG_LOOKUP, data=payload, headers=self.header)
        else:
            request = self.session.post(IMG_DATA, data=payload, headers=self.header)
        if self._proceed_with_download(local_path, request=request):
            with open(local_path, 'w') as OUT:
                OUT.write(request.text)
            LOG.info("Downloading {} for {} to\n    {}\n\n".format(datatype, prefix, local_path))
            return "successful"
        else:
            return "skipped"

    def _download_data(self, url, local_path):
        """ 
        Downloads a file in chunks.
        Returns: Nothing
        Excepts: Nothing
        """

        #download the file in chunks
        LOG.info("Downloading {} to\n    {}\n\n".format(url, local_path))
        request = self.session.get(url, stream=True)
            
        try:
            with open(local_path, 'wb') as fh:
                for chunk in request.iter_content(chunk_size=1024):
                    if chunk:  # filter out keep-alive new chunks
                        fh.write(chunk)
                        fh.flush()
        except IOError:

            LOG.warning("Problem downloading to local path {}".format(local_path))

    def _proceed_with_download(self, local_path, file=None, request=None):
        """ 
        Checks the local path to ensure no unintentional clobbering.
        
        Note, clobbering can still occur between the time of check and the time of 
        write if there are unknown race conditions.
        """

        LOG.debug("Checking if we should proceed with download...")
        if os.path.exists(local_path):
            if self.force_overwrite:
                LOG.warning("Force overwriting: {}".format(local_path))
                self.overwritten_paths.append(local_path)
                return True
            else:
                if self.resume:
                    #### TODO implement file checking here, first MD5, then filesize
                    # for IMG, check size?
                    LOG.info("Found file: {}. Skipping Download!".format(local_path))
                    self.skipped_paths.append(local_path)
                    return False
                else:
                    
                    LOG.warning("Neither force overwrite nor resume are enabled. Refusing to download to {}".format(local_path))
                    return False

        else:
            return True




def main2(args):
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

    """
    Sets up the Interface and makes an Organism for each id and returns a list of Organisms.
    """
   
    # set up the interface

def init_organisms(interface, ids, id_type, names = None):
    """ Returns a list of organisms for a list of ids and optional names. """

    organisms = []
    for tmp_id in ids:
        if id_type == "taxon_oid":
            org = JGIOrganism(interface, taxon_oid=tmp_id)
        elif args.id_type == "proj_id":
            org = JGIOrganism(interface, proj_id=tmp_id)
        organisms.append(org)

    # add names if present
    if names:
        for org, name in zip(organisms, names):
            org.prefix = name

    return organisms

def download_data_for_all_organisms(organisms, datatype):
    """
    Downloads a piece of data for all organisms.
    Powers through all exceptions.
    """
    for organism in organisms:
        try:
            organism.download_data(datatype)
        except (AttributeLookupError, PortalError) as e:
            LOG.error("Organism {} failed datatype '{}'\n".format(str(organism), datatype) + e.args[0])
            continue

def standard_pipeline(args):

    interface = JGIInterface(login_file=args.login, force_overwrite=args.force, resume=args.resume, newest_only=args.newest)

    if args.download:
        
        local_path = args.download.split("/")[-1]
        interface._download_data(args.download, local_path)
        sys.exit()
    
    else:
        ids = _read_file_or_split_arg(args.ids)
        if args.names:
            names = _read_file_or_split_arg(args.names)
        else:
            names = None
        organisms = init_organisms(interface, ids, args.id_type, names)

        # if the manual option has been passed, exit here
        if args.manual:
            return organisms

        # check for args.get
        if args.get:
            for datatype in args.get:
                download_data_for_all_organisms(organisms, datatype)

        # if not args.get we want to list available files
        else:
            for organism in organisms:
                organism.print_available_files(portal="both")

        print("PRINTING DOWNLOAD REPORT")
        for organism in organisms:
            print(str(organism))
            print("--------------------")
            organism.print_download_report()
            print(end="\n\n")

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
        [items.append(item) for item in arg.split(delim)]

        
    return items

def arg_type_get(arg):
    """
    Turns a text string input from the command line to a list of touples like the program requires.
    """

    # see if get is for IMG data, will check for id type later
    if arg in JGIOrganism.IMG_DATA:
        return arg

    try:
        (folder, regex) = arg.split("|", 1)
        return folder, regex
    except:
        raise
        #raise argparse.ArgumentTypeError("Error processing the --get argument.")

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

    """, choices=['proj_id', 'taxon_oid'])
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
    parser.add_argument("-d", "--download", help="""
single path (as URL) to download

    """)
    parser.add_argument("-f", "--force", help="""
forces download of files, overwrites files that already have a local path

otherwise, program will not overwrite paths.

    """, action="store_true")
    parser.add_argument("-n", "--newest", help= """
set to keep the newest file only when multiple files match the same regex
    """, action="store_true")
    parser.add_argument("-r", "--resume", help="""
skips downloading if local file is found

if this option is not given and a duplicate file is encountered
without the --force option, program will exit

file must have the same name to be considered whole file.

eventually plan to implement size or MD5 checks to ensure
the file is complete before skipping

    """, action="store_true")
    parser.add_argument("-v", help="""
sets the verbosity level for the logging module

    """, choices=["DEBUG", "INFO", "WARNING"], default="INFO")
    parser.add_argument("-m", "--manual", help="""
don't use the automated workflow, just give me a a variable
"organisms" with the organisms from my ids. Meant to be used
for interactive python usage
    """, action="store_true")
    
   
    args = parser.parse_args()

    # set the logging level
    logging.basicConfig(level=args.v)

    organisms = standard_pipeline(args)
