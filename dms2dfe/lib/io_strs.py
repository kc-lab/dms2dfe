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
