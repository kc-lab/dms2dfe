#!usr/bin/python

# Copyright 2016, Rohan Dandage <rraadd_8@hotmail.com,rohan@igib.in>
# This program is distributed under General Public License v. 3.  

"""
================================
``io_stats``
================================
"""
import pandas as pd
import numpy as np
import logging
from statsmodels.stats.weightstats import DescrStatsW,CompareMeans
from statsmodels.sandbox.stats.multicomp import multipletests

from dms2dfe.lib.io_dfs import debad 

def testcomparison(df,smp1_cols,smp2_cols,test='ttest'):
    col_stat='stat %s' % test
    col_pval='pval %s' % test
    df.loc[:,col_stat]=np.nan
    df.loc[:,col_pval]=np.nan
    for i in df.index:
        X=DescrStatsW(df.loc[i,smp1_cols].as_matrix())
        Y=DescrStatsW(df.loc[i,smp2_cols].as_matrix())
        if test=='ttest':
            df.loc[i,col_stat],df.loc[i,col_pval],tmp=CompareMeans(X,Y).ttest_ind()
        if test=='ztest':
            df.loc[i,col_stat],df.loc[i,col_pval]=CompareMeans(X,Y).ztest_ind()
    return df

from scipy import stats
from dms2dfe.lib.io_dfs import denan
def get_r2(data,xcol,ycol,log=None):
    data=denan(data.loc[:,[xcol,ycol]],axis='rows',condi='any')
    if len(data)!=0:
        if not log is None:
            if log==2:
                data=debad(data,axis=0,condi='any',bad=0)                
                data=np.log2(data)
                data=debad(data,axis=0,condi='any',bad='nan')
        slope, intercept, r_value, p_value, std_err = stats.linregress(data.loc[:,xcol],data.loc[:,ycol])
        return r_value
    else:
        logging.error("one/both cols are empty")
        return 0