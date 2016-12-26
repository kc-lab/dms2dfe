#!usr/bin/python

# Copyright 2016, Rohan Dandage <rraadd_8@hotmail.com,rohan@igib.in>
# This program is distributed under General Public License v. 3.  

"""
================================
``io_dfs``
================================
"""

import pandas as pd

def set_index(data,col_index):
    if col_index in data:
        data=data.reset_index().set_index(col_index)
        if 'index' in data:
            del data['index']
        return data
    elif data.index.name==col_index:
        return data


# dfs
def concat_cols(df1,df2,idx_col,df1_cols,df2_cols,
                df1_suffix,df2_suffix):
#     df1_cols=[ for col in df1_cols]
    df1=df1.set_index(idx_col)
    df2=df2.set_index(idx_col)    
    combo=pd.concat([df1.loc[:,df1_cols],df2.loc[:,df2_cols]],axis=1)
    # find common columns and rename them
    common_cols=[col for col in df1_cols if col in df2_cols]
    for col in common_cols:
        df1_cols[df1_cols.index(col)]="%s%s" % (col,df1_suffix)
        df2_cols[df2_cols.index(col)]="%s%s" % (col,df2_suffix)
    combo.columns=df1_cols+df2_cols
    combo.index.name=idx_col
    return combo

def get_colmin(data):
    data=data.T
    colmins=[]
    for col in data:
        colmins.append(data[col].idxmin())
    return colmins
