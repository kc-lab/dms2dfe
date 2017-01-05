#!usr/bin/python

# Copyright 2016, Rohan Dandage <rraadd_8@hotmail.com,rohan@igib.in>
# This program is distributed under General Public License v. 3.  

"""
================================
``io_strs``
================================
"""


def isstrallowed(s,form):
    """
    Checks is input string conforms to input regex (`form`).

    :param s: input string.
    :param form: eg. for hdf5: `"^[a-zA-Z_][a-zA-Z0-9_]*$"`
    """
    import re
    match = re.match(form,s)
    return match is not None

def convertstr2format(col,form):
    """
    Convert input string to input regex (`form`).
    
    :param col: input string.
    :param form: eg. for hdf5: `"^[a-zA-Z_][a-zA-Z0-9_]*$"`
    """
    if not isstrallowed(col,form):
        col=col.replace(" ","_") 
        if not isstrallowed(col,form):
            chars_disallowed=[char for char in col if not isstrallowed(char,form)]
            for char in chars_disallowed:
                col=col.replace(char,"_")
    return col

def make_pathable_string(s):
    return s.replace(" ","_").replace("$","").replace("\\","")\
            .replace("(","").replace(")","")\
            .replace("{","").replace("}","").replace("-","")\
            .replace(":","").replace("^","").replace("+","").replace("'","").replace("\"","")\
            .replace("\n","_").replace("\t","_")

def get_logger(argv=None):
    import logging
    import datetime
    log_format='[%(asctime)s] %(levelname)s\tfrom %(filename)s in %(funcName)s(..):%(lineno)d: %(message)s'
    logging.basicConfig(format=log_format,
                        level=logging.DEBUG,)
    if not argv is None:
        log_fh="%s_%s" % (make_pathable_string(str(datetime.datetime.now())),'_'.join(argv).replace('/','|'))
        print log_fh
        logging.basicConfig(filename=log_fh)
        console = logging.StreamHandler()
        console.setLevel(logging.INFO)
        formatter = logging.Formatter(log_format)
        console.setFormatter(formatter)
        logging.getLogger('').addHandler(console)
    # logging.info('#START')
    return logging
