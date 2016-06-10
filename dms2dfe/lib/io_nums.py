#!usr/bin/python

# Copyright 2016, Rohan Dandage <rraadd_8@hotmail.com,rohan@igib.in>
# This program is distributed under General Public License v. 3.  

"""
================================
``io_nums``
================================
"""

def is_numeric(obj):
    """
    This detects whether an input object is numeric or not.
    """
    try:
        obj+obj, obj-obj, obj*obj, obj**obj, obj/obj
    except ZeroDivisionError:
        return True
    except Exception:
        return False
    else:
        return True
    
def str2num(x):
    """
    This extracts numbers from strings. eg. 114 from M114R.
    """
    return int(''.join(ele for ele in x if ele.isdigit()))