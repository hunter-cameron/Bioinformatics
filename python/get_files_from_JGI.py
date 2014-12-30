#!/usr/bin/python

'''
This script allows files to be downloaded from JGI in an automated mannor. It uses
k-mer text matching to align filenames because JGI often does not keep
consistent naming conventions
'''

import requests     #not in core python. install using 'pip install --user requests'
import warnings
import re
import os
import sys
import getopt
import time

global LOGIN, XML, DATA, LOOKUP
LOGIN = 'https://signon.jgi.doe.gov/signon/create'
XML = 'http://genome.jgi.doe.gov/ext-api/downloads/get-directory?organism='     #complete with organism name
DATA = 'http://genome.jgi.doe.gov/'     #complete with path to file from the organism XML
LOOKUP = "http://genome.jgi.doe.gov/lookup?"


def usage():
    print ""
    print "Usage: automate_jgi.py -l <file> -i <file> -o <file> -g 'JGI|raw' (--login, --ids, --out, --get)"
    print "         -l, --login <file>          File with login credentials.(username\\npassword)"
    print "         -i, --ids <file>            File with a list of JGI Taxon Ids"
    print "         -o, --out <file>            Prefix path for output"
    print "         -g, --get <'JGI|raw'>       Folder to download. Currently supports either JGI data or raw data."
    print ""


class JGI_genome:

    def __init__(self, parent, proj_id):
        self.proj_id = proj_id
        self.parent = parent
        #print parent

        self._get_name()
        self._parse_XML()

    def _get_name(self):

        my_params = {'keyName': 'jgiProjectId', 'keyValue': self.proj_id}
        name_html = self.parent.session.get(LOOKUP, params=my_params)    #use the parent session to maintain cookies
        #print name_html.url
        match = re.search('href="/(.*)/.*\.info\.html', name_html.text)     #search for href="/(NAME)/ANYTHING.info.html
        #print match.group(1)
        #check for a match
        if match:
            self.name = match.group(1)       #get the first matched group
            self._get_XML()
        else:
            sys.exit("Name not found for id: %s" % self.proj_id)

            return
            #instead of return, might be nice to remove the reference in the parent list


    def _get_XML(self):

        link = XML + self.name
        #print "XML:", link 
        self.xml = self.parent.session.get(link).text
        #print self.xml
        self._parse_XML()

    #parse out the links in a 2 step approach
    #1. get the folder
    #2. parse out all links in the folder
    def _parse_XML(self):



        #match = re.findall('<folder\ name="IMG\ Data">(.*?)</folder>', self.xml, re.DOTALL)
        match = re.search('<folder\ name="IMG\ Data">(.*?)</folder>', self.xml, re.DOTALL)
        if match:
            #print "MATCH: IMG", match
            #print match.group(1)
            self.IMG_data = re.findall('url=(\/.*?\.tar\.gz)"', match.group(1))     #original regex had a " after the =; but this returns an invalid path (the whole one JGI gives)
            #print self.IMG_data[1]
            #sys.exit()
        else:
            #do something
            warnings.warn("IMG Link not found for %s" % self.proj_id)
            sys.exit(2)

        match = None
        #print self.xml
        match = re.search('<folder\ name="Raw\ Data">(.*?)</folder>', self.xml, re.DOTALL)
        #print match.group(1)
        if match:
            print match.group(1)
            self.raw_data = re.findall('url=(\/.*?\.gz)"', match.group(1))     #original regex had a " after the =; but this returns an invalid path
        else:
            #do something else
            #warnings.warn("Link not found for %s" % self.proj_id)
            sys.exit("Raw link not found for %s" % self.proj_id)

    def download_data(self, dtype, prefix):
        #mapping function to make each path absolute
        #def append_url(x): return DATA + x

        print "JGI", self.IMG_data
        print "raw", self.raw_data



        if (dtype == "IMG"):
            urls = self.IMG_data
        elif (dtype == "raw"):
            urls = self.raw_data
        else:
            sys.exit("Incorrect datatype")

        #download all the urls
        for url in urls:
            url = DATA + url
            print url
            local_filename = prefix + "/" + url.split('/')[-1]

            #download the file in chunks
            r = self.parent.session.get(url, stream=True)
            with open(local_filename, 'wb') as fh:
                for chunk in r.iter_content(chunk_size=1024):
                    if chunk: # filter out keep-alive new chunks
                        fh.write(chunk)
                        fh.flush()


class JGI_interface:


    def __init__(self, login_file, jgi_file):

        self.login(login_file)
        self._make_genomes(jgi_file)


    def login(self, login_file):

        #parse the login info
        if login_file == '':
            print "Login file not specified"
            usage()
            sys.exit(2)

        with open(login_file, 'r') as file_handle:
            login_data = {
                'login': file_handle.readline().rstrip(),
                'password':  file_handle.readline().rstrip()
            }

        self._authenticate(login_data)

    def get_IMG_data(self, prefix):
        x = 1
        for genome in self.JGI_genomes:
            print "Downloading Genome: %s" % x
            x += 1
            genome.download_data("IMG", prefix)

    def get_raw_data(self, prefix):
        x = 1
        for genome in self.JGI_genomes:
            print "Downloading Genome: %s of %s." % (x, len(self.JGI_genomes))
            genome.download_data("raw", prefix)


######Internal methods

    def _authenticate(self, login_data):

        #make a session to store cookie information like a browser would
        self.session = requests.session()

        #login
        r = self.session.post(LOGIN, data=login_data)

        match = re.search('You have signed in successfully.', r.text)
        if (not match):
            print "Login may not have been successful"
            time.sleep(2)
        #print r.text       #uncomment to display raw login info. Useful to look through if having problems 

    def _make_genomes(self, jgi_file):

        self.JGI_genomes = []

        #read the jgi_file for proj_ids
        with open(jgi_file, 'r') as fh:
            for line in fh:
                line = line.rstrip()
                fields = line.split("\t")
                proj_id = fields[0]
                proj_id = JGI_genome(self, proj_id)
                self.JGI_genomes.append(proj_id)




'''
BEGIN MAIN USAGE
'''



if (__name__ == "__main__"):

    os.nice(20)     #give this program lowest priority to give up CPU to others while downloads are going on
    login_file = ''
    ids_file = ''
    dfile = ""
    prefix = os.getcwd()

    #getopt returns tuples - but what is stored into the args variable?
    try:
        opts, args = getopt.getopt(sys.argv[1:], "l:i:o:g:", ["login=","ids=", "out", "get="])
    except getopt.GetoptError as err:
        print str(err)
        usage()
        sys.exit(2)

    #looks like black magic - opts contains tuples which are a list of values (eg. "file", "arg") like database entries that stay together
    #this for loop splits opts into the respective tuples
    for opt, arg in opts:
        if opt in ("-l", "--login"):
            login_file = arg
        elif opt in ("-i", "--ids"):
            ids_file = arg
        elif opt in ("-g", "--get"):
            if (arg == "JGI" or arg == "raw"):
                dfile = arg
            else:
                usage()
                sys.exit("Inappropriate value: %s for -g, --get" % arg)
        elif opt in ("-o", "--out"):
            prefix = arg



    interface = JGI_interface(login_file, ids_file)


    if (dfile == "JGI"):
        interface.get_IMG_data(prefix)
    elif (dfile == "raw"):
        interface.get_raw_data(prefix)



