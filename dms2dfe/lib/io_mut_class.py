# Copyright 2017, Rohan Dandage <rraadd_8@hotmail.com,rohan@igib.in>
# This program is distributed under General Public License v. 3.  

"""
================================
``io_mut_files``
================================
"""
from __future__ import division
import sys
from os.path import splitext,exists,basename,abspath,dirname
from os import makedirs,stat
import pandas as pd
import numpy as np
from glob import glob
import logging
import subprocess

def get_2SD_cutoffs(d,reps,N=False):
    t=d.copy()
    t.loc[:,'reps']=t.loc[:,reps[0]]-t.loc[:,reps[1]]
    if N and ('mut' in t):
        t.loc[(t.loc[:,'mut']==t.loc[:,'ref']),'reps']=np.nan
    mu=t.loc[:,'reps'].mean()
    sigma=t.loc[:,'reps'].std()
    return mu+sigma*2,mu,sigma

def get_repli_FCN(d,csel='.NiA_tran.sel',cref='.NiA_tran.ref'):
    sels=np.sort([c for c in d.columns if csel in c])
    refs=np.sort([c for c in d.columns if cref in c])
#     sel1=[c for c in sels if 'replicate_1' in c][0]
#     sel2=[c for c in sels if 'replicate_2' in c][0]
#     ref1=[c for c in refs if 'replicate_1' in c][0]
#     ref2=[c for c in refs if 'replicate_2' in c][0]
    
#     for rep in ['replicate_1','replicate_2']:
    FCAi=1
    cols_FCA=[]
    for refi,ref in enumerate(refs):
        for seli,sel in enumerate(sels):
            if refi==seli:
                coln='FCA%s_reps' % FCAi
                d.loc[:,coln]=d.loc[:,sel]-d.loc[:,ref]
                d.loc[(d.loc[:,'mut']==d.loc[:,'ref']),coln]=np.nan
                cols_FCA.append(coln)
                FCAi+=1
#     FCAi=1
#     cols_FCA=[]
#     for refi,ref in enumerate(refs):
#         for seli,sel in enumerate(sels):
#             if refi>=seli:
#                 coln='FCA%s' % FCAi
#                 d.loc[:,coln]=d.loc[:,sel]-d.loc[:,ref]
#                 d.loc[(d.loc[:,'mut']==d.loc[:,'ref']),coln]=np.nan
#                 cols_FCA.append(coln)
#                 FCAi+=1
    return d.loc[:,cols_FCA]

def data_fit2cutoffs(d,sA,sB,N=True):
    refs=[c for c in d.columns if sA in c]
    sels=[c for c in d.columns if sB in c]
    _,mu1,sigma1=get_2SD_cutoffs(d,refs)
    _,mu2,sigma2=get_2SD_cutoffs(d,sels)
    return np.mean([mu1,mu2])+2*np.sqrt(np.mean([sigma1,sigma2]))

def class_fit(d,col_fit='FiA',FC=True,zscore=False): #column of the data_fit
    """
    This classifies the fitness of mutants into beneficial, neutral or, deleterious.
    
    :param d: dataframe of `data_fit`.
    :returns d: classes of fitness written in 'class-fit' column based on values in column 'FiA'. 
    """
    cols_reps=[c for c in dctrl if '.NiA_tran.ref' in c]
    if FC and (len(cols_reps)==2):
        up,_,_=get_2SD_cutoffs(d,cols_reps,N=True)
        dw=-1*up
    else:
        if zscore:
            up,dw=-2,2
        else:
            up,dw=0,0
    d.loc[d.loc[:,col_fit]>+up,    'class_fit']="enriched"
    d.loc[((d.loc[:,col_fit]>=dw) & (d.loc[:,'FiA']<=up)),'class_fit']="neutral"
    d.loc[d.loc[:,col_fit]<dw,    'class_fit']="depleted"
    return d

def class_comparison(dA,dB):
    """
    This classifies differences in fitness i.e. relative fitness into positive, negative or robust categories. 
    
    :param dc: dataframe with `dc`. 
    :returns dc: dataframe with `class__comparison` added according to fitness levels in input and selected samples in `dc`
    """
    dc=get_repli_FCN(dA).join(get_repli_FCN(dB),lsuffix='_ctrl',rsuffix='_test')
    up=data_fit2cutoffs(dc,sA='_reps_test',sB='_reps_ctrl',N=False)
    dw=-1*up

    diff=dB.loc[:,'FCN']-dA.loc[:,'FCN']
    diff=diff.reset_index()
    mutids_up=diff.loc[(diff.loc[:,'FCN']>up),'mutids'].tolist()
    mutids_dw=diff.loc[(diff.loc[:,'FCN']<dw),'mutids'].tolist()

    d.loc[d.loc[:,col_fit]>+up,    'class_comparison']="enriched"
    d.loc[((d.loc[:,col_fit]>=dw) & (d.loc[:,'FiA']<=up)),'class_comparison']="neutral"
    d.loc[d.loc[:,col_fit]<dw,    'class_comparison']="depleted"
    return dc
