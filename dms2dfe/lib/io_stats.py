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
from statsmodels.stats.weightstats import DescrStatsW,CompareMeans
from statsmodels.sandbox.stats.multicomp import multipletests

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
