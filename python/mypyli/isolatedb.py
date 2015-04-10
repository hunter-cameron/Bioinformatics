
from __future__ import print_function
import pandas

# MAGIC
ISOLATE_DB = '/nas02/home/h/j/hjcamero/scripts/python/mypyli/isolate_db'



def read_database(db_f):
    """
    make a database of style: [[col1], [col2], [col3]] where the header is the first indedx
    in each col
    """
    with open(db_f, 'r') as IN:
        headers = IN.readline()
        database = [[name] for name in headers[:-1].split("\t")]
        
        for row in IN:
            for indx, value in enumerate(row[:-1].split("\t")):
                database[indx].append(value)

    return database

database = read_database(ISOLATE_DB)


def convert_values(values, conv_from='taxon_id', conv_to='isolate'):
    global database
    #print(database)
    # I'll need to wrap each of these with try except 
    from_indx = [col[0] for col in database].index(conv_from)
    to_indx = [col[0] for col in database].index(conv_to)

    conv_values = []
    for value in values:
        value = str(value)
        val_indx = database[from_indx].index(value)

        conv_values.append(database[to_indx][val_indx])

    return conv_values

def convert_dataframe(frame):
    indices = frame.index.values
    new_indx = convert_values(indices, conv_from='taxon_id', conv_to='isolate')
    #print(new_indx)
    new_frame = frame.copy(deep=True)
    new_frame['index'] = new_indx

    new_frame.set_index('index', inplace=True)

    #print(new_frame)
    return new_frame

