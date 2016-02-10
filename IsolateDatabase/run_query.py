
import argparse
import getpass
import pymysql
import os
import sys

import logging
logging.basicConfig()
LOG = logging.getLogger(__name__)
LOG.setLevel("INFO")

def connect_to_database(host, user, password, db):
    """ Connects to a database and returns a connection. """
    LOG.info("Attempting to connect to database...")
    conn = pymysql.connect(host=host, user=user, passwd=password, db=db)
    LOG.info("Connected!")

    return conn

def get_sql_statement(sql):
    """ Gets a sql statement from either a statement or a file """
    if os.path.isfile(sql):
        with open(sql, 'r') as IN:
            stmt = IN.read()
    else:
        stmt = str(sql)

    stmt = stmt.split(";\n")

    return stmt

def run_sql_statement(conn, stmt, header=False):
    """ Runs a sql statement and iterates over the output """
    LOG.info("Running sql query...")
    cur = conn.cursor()
    
    # only return results for the last executed statement 
    for cmd in stmt:
        LOG.info("Running: '{}'".format(cmd))
        cur.execute(cmd)

    LOG.info("Yielding results.")

    if header is True:
        yield "\t".join([str(field[0]) for field in cur.description])

    for row in cur.fetchall():
        yield "\t".join([str(field) for field in row])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Connects to a MySql database and runs a query and dumps the results to stdout")
    parser.add_argument("-host", help="host to connect to [%(default)s]", default="127.0.0.1")
    parser.add_argument("-user", help="the MySql user")
    parser.add_argument("-login", help="file to parse the username and login from")
    parser.add_argument("-db", help="the name of the database to connct to")
    parser.add_argument("-sql", help="SQL statement or file to execute", required=True)
    parser.add_argument("-header", help="print out a header first [%(default)s]", action='store_true')
    args = parser.parse_args()
    

    if args.login:
        with open(args.login, 'r') as IN:
            args.user = IN.readline().strip()
            password = IN.readline().strip()
    else:
        if args.user:
            password = getpass.getpass("Password for '{}': ".format(args.user))
        else:
            raise ValueError("No username or login file supplied!")

    conn = connect_to_database(host=args.host, user=args.user, password = password, db=args.db)
    stmt = get_sql_statement(args.sql)

    for row in run_sql_statement(conn, stmt, args.header):
        print(row)

