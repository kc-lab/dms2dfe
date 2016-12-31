#!usr/bin/python

# Copyright 2016, Rohan Dandage <rraadd_8@hotmail.com,rohan@igib.in>
# This program is distributed under General Public License v. 3.  

"""
================================
``io_mut_files``
================================
"""
from __future__ import division
import sys
from os.path import splitext,exists,basename
from os import makedirs,stat
import pandas as pd
import numpy as np
from glob import glob
import logging
logging.basicConfig(format='[%(asctime)s] %(levelname)s\tfrom %(filename)s in %(funcName)s(..):%(lineno)d: %(message)s',level=logging.DEBUG) # filename=cfg_xls_fh+'.log'
from dms2dfe.lib.fit_curve import fit_gauss_params
from dms2dfe.lib.global_vars import aas_21,cds_64,mut_types_form,mut_types_NorS
from dms2dfe.lib.convert_seq import cds2aas
from dms2dfe.lib.io_seq_files import getdepth_cds,getdepth_ref
from dms2dfe.lib.io_data_files import getusable_lbls_list,getusable_fits_list,getusable_comparison_list
from dms2dfe.lib.io_dfs import concat_cols
from dms2dfe.lib.io_nums import str2num

def makemutids(data_lbl,refis):
    """
    This makes the mutation ids eg. M001T or A024T from mutation matrix.
    
    :param data_lbl: dataframe with columns of reference (`ref`) and mutation (`mut`).
    :param refis: index of reference amino acids/ codons.
    :returns mutids: list of mutation ids.
    """
    mutids=[]
    data_lbl=data_lbl.reset_index()
    reflen=len(refis)
    for muti in range(len(data_lbl)//reflen):
        refii=0
        for refi in refis:
            mutids.append("%s%03d%s" % (data_lbl.loc[muti*reflen+refii,'ref'],refi,data_lbl.loc[muti*reflen+refii,'mut']))
            refii+=1
    if len(data_lbl)!=len(mutids):
        logging.error("len(data_lbl)!=len(mutids) bcz %d != %d" % (len(data_lbl),len(mutids)))
        sys.exit()
    return mutids    

def makemutids_fromprtseq(prt_seq,muts=None):
    refs=list(prt_seq)
    refis=range(1,len(refs)+1,1)
    if muts is None:
        from dms2dfe.lib.global_vars import aas_21
        muts=aas_21
    mutids=[]
    for i in range(len(refs)):
        for mut in muts: 
            mutid="%s%03d%s" % (refs[i],refis[i],mut)
            mutids.append(mutid.replace('*','X'))
    return mutids


def mutids_converter(mutids,out,type_form):
    if type_form=='aas':
        offset=0
    elif type_form=='cds':
        offset=2
        
    if out=='refi':
        return [str2num(mutid) for mutid in mutids]    
    elif out=='ref':
        return [mutid[:1+offset] for mutid in mutids]
    elif out=='mut':
        return [mutid[-1-offset:] for mutid in mutids]
    elif out=='refrefi':
        return [mutid[:4+offset] for mutid in mutids]

def mat_cds2mat_aas(mat_cds,host) :
    """
    This converts mutation matrix of codon to that of amino acid.
    
    :param mat_cds: codon level mutation matrix.
    :param host: name of host organism for choosing codon table. [coli | yeast | sapiens].
    :returns mat_aas: amino acid level mutation matrix.
    """
    aas_wt=[cds2aas(cd,host) for cd in list(mat_cds.index)]
    aas_64=[cds2aas(cd,host) for cd in list(mat_cds.columns)]
    mat_aas_64 =mat_cds.fillna(0)
    mat_aas_64.columns=aas_64
    mat_aas=mat_aas_64.groupby(mat_aas_64.columns, axis=1).sum()
    mat_aas    =mat_aas.loc[:,aas_21]
    mat_aas.index=aas_wt
    mat_aas.columns.name='mut'
    mat_aas.index.name='ref'
    return mat_aas

def getNS(data_all) : # data (col): can be cds or aas
    """
    This function separates the non-synonymous(N) and synonymous(S) mutations.
    
    :param data_all: dataframe with frequency of mutations (`data_lbl`).
    :returns data_NS: dataframe with non-synonymous(N) and synonymous(S) mutations.
    """
    data_NS=pd.DataFrame(columns=['NiS','NiN'],index=data_all.index);#data_NS.loc[:,'mutids']=data_all.loc[:,'mutids']
    #data_S=pd.DataFrame(columns=data_all.columns,index=data_all.index);data_S.loc[:,'mutids']=data_all.loc[:,'mutids']
    rowcount=0
    for rowi,row in data_all.iterrows() :
        if rowi[0]== rowi[1] :
            data_NS.ix[rowcount,'NiS']=data_all.ix[rowcount,'NiA']    
        else :
            data_NS.ix[rowcount,'NiN']=data_all.ix[rowcount,'NiA']
        rowcount+=1    
    return data_NS

def collate_cctmr(lbl_mat_cds,cctmr):    
    """
    If the reference sequence is a concatamer, this function collates the repeating mutation matrix. 
    
    :param lbl_mat_cds: codon level mutation matrix.
    :param cctmr: tuple with (1-based) boundaries of concatamers.
    :returns lbl_mat_cds: collated codon level mutation matrix.
    """
    if (lbl_mat_cds.index.values[(cctmr[0][0]-1)] == lbl_mat_cds.index.values[(cctmr[1][0]-1)]) and \
    (lbl_mat_cds.index.values[(cctmr[0][1]-1)] == lbl_mat_cds.index.values[(cctmr[1][1]-1)]):
        lbl_mat_cds_cctmr1=lbl_mat_cds.iloc[(cctmr[0][0]-1):(cctmr[0][1]-1),:]
        lbl_mat_cds_cctmr2=lbl_mat_cds.iloc[(cctmr[1][0]-1):(cctmr[1][1]-1),:]
        lbl_mat_cds=lbl_mat_cds_cctmr1.fillna(0)+lbl_mat_cds_cctmr2.fillna(0)
        lbl_mat_cds.loc[:,"refi"]=lbl_mat_cds_cctmr1.loc[:,"refi"]
    else:
        logging.error("cctmr do not conform %s!=%s or %s!=%s" % (lbl_mat_cds.index.values[(cctmr[0][0]-1)], \
                                                                 lbl_mat_cds.index.values[(cctmr[1][0]-1)], \
                                                                 lbl_mat_cds.index.values[(cctmr[0][1]-1)], \
                                                                 lbl_mat_cds.index.values[(cctmr[1][1]-1)]))
    return lbl_mat_cds.replace(0,np.nan)

def repli2avg(replis,prj_dh,type_form):
    """
    This averages the mutation data of replicates (data_lbl format).
    
    :param replis: list of names of replicate samples.
    :param prj_dh: path to project directory.
    :param type_form: type of mutations codon (`cds`) or amino acid (`aas`) form.
    :returns data_avg: dataframe with averaged values.
    """
    data_avg_ar=np.array([])
    for replii in replis :
        data_replii=pd.read_csv("%s/data_lbl/%s/%s" % (prj_dh,type_form,replii))
        if not 'mut' in data_replii:
            data_replii.loc[:,'mut']\
            =mutids_converter(data_replii.loc[:,'mutids'],
                              'mut',type_form)
        if not 'ref' in data_replii:
            data_replii.loc[:,'ref']\
            =mutids_converter(data_replii.loc[:,'mutids'],
                              'ref',type_form)
        
        data_replii=data_replii.set_index(['mut','ref','mutids'])
        data_replii=data_replii.fillna(0)
        if len(data_avg_ar)==0 :
            data_avg_index  =data_replii.index.values #fill with zeros                        
            data_avg_index = [  np.array([t[0] for t in data_avg_index ]) ,\
                                np.array([t[1] for t in data_avg_index ]) ,\
                                np.array([t[2] for t in data_avg_index ])]
            data_avg_columns=data_replii.columns.values
            data_avg_ar=data_replii.as_matrix()
        else :
            data_avg_ar=(data_avg_ar+data_replii.as_matrix()) /2
    data_avg=pd.DataFrame(data_avg_ar, index=data_avg_index, columns=data_avg_columns)
    data_avg.index.names=['mut','ref','mutids']
    data_avg=data_avg.reset_index()
    data_avg=data_avg.replace(0, np.nan)
    return data_avg
                
def repli2data_lbl_avg(prj_dh):    
    """
    This averages frequency data of replicates in data_lbl.
    
    :param prj_dh: path to project directory.
    """
    from dms2dfe.lib.global_vars import mut_types_form

    if exists('%s/cfg/repli' % prj_dh):
        repli=pd.read_csv('%s/cfg/repli' % prj_dh)
        repli=repli.set_index('varname')
        for avg_lbl,replis in repli.iterrows() :
            if not pd.isnull(avg_lbl):
                replis=replis[~replis.isnull()]
                replis=list(replis)
                for type_form in mut_types_form : # get aas or cds 
                    if not exists("%s/data_lbl/%s/%s" % (prj_dh,type_form,avg_lbl)):
                        repli_keys=[("data_lbl/%s/%s" % (type_form,replii)) for replii in replis]                 
                        if all(exists("%s/%s" % (prj_dh,repli_key)) for repli_key in repli_keys):
                            logging.info("processing : %s/%s" % (type_form,avg_lbl))
                            data_avg=repli2avg(replis,prj_dh,type_form) # list
                            data_avg.reset_index().to_csv('%s/data_lbl/%s/%s' % (prj_dh,type_form,avg_lbl),index=False)  
                        else: 
                            for repli_key in repli_keys:
                                if not exists("%s/%s" % (prj_dh,repli_key)):
                                    logging.warning("can not find: %s/%s" % (prj_dh,repli_key))
                    # else:
                    #     logging.info("already processed %s" % avg_lbl)
    else:
        logging.warning("skipping repli2data_lbl_avg")

def mut_mat_cds2data_lbl(lbli,lbl_mat_mut_cds_fh,
                    host,prj_dh,reflen,cctmr,Ni_cutoff,fsta_fh,clips=None):
    """
    This function converts mutation matrix (produced from .mat) file into data_lbl format.
    
    :param lbli: name of the sample.
    :param lbl_mat_mut_cds_fh: path to codon level mutation matrix (.mat_mut_cds).
    :param host : name of host organism for choosing codon table. (coli | yeast | sapiens).
    :param prj_dh: path to project directory.
    :param reflen: nucleotide level length of the reference sequence.
    :param cctmr: tuple with (1-based) boundaries of concatamers.
    :param Ni_cutoff: threshold of depth per codon to be processed further. 
    """
    if stat(lbl_mat_mut_cds_fh).st_size != 0 :
        lbl_mat_cds=pd.read_csv(lbl_mat_mut_cds_fh) # get mat to df
        lbl_mat_cds=lbl_mat_cds.set_index("ref_cd",drop=True) # make codons as the index
        if 'Unnamed: 0' in lbl_mat_cds.columns:
            del lbl_mat_cds['Unnamed: 0']
        lbl_mat_cds.columns.name='mut'
        lbl_mat_cds.index.name='ref'
        if cctmr != None:
            reflen=(cctmr[0][1]-1)*3
            lbl_mat_cds=collate_cctmr(lbl_mat_cds,cctmr)
        lbl_mat_aas=mat_cds2mat_aas(lbl_mat_cds,host)# convert to aa data
        # lbl_mat_cds.to_csv('test')
        for type_form in mut_types_form : # get aas or cds  
            if not exists("%s/data_lbl/%s/%s" % (prj_dh,type_form,str(lbli))):
                # print lbl_mat_aas.tail(10)
                if type_form=="aas":
                    data_lbl=pd.DataFrame(lbl_mat_aas.unstack())
                elif type_form=="cds":
                    data_lbl=pd.DataFrame(lbl_mat_cds.drop("refi",axis=1).unstack())
                # print data_lbl.iloc[170:180,:]
                ref_prt_len=(float(reflen)//3)
                if (len(data_lbl)//ref_prt_len==21) or (len(data_lbl)//ref_prt_len==64):
                    data_lbl.columns=["NiA"]  # rename col Ni  
                    data_lbl.loc[:,'mutids']=makemutids(data_lbl,lbl_mat_cds.loc[:,"refi"]) # add mutids
                    data_lbl=pd.concat([data_lbl,getNS(data_lbl)], axis=1)
                    
                    for type_NorS in mut_types_NorS : # N OR S  
                        data_lbl.loc[:,('Ni%scut' % type_NorS)]=data_lbl.loc[:,'Ni%s' % type_NorS]
                        data_lbl.loc[data_lbl.loc[:,('Ni%scut' % type_NorS)]<Ni_cutoff,('Ni%scut' % type_NorS)]=np.nan \
                        #Ni_cutoff=8 varscan default
                        data_lbl.loc[:,('Ni%scutlog' % type_NorS)]\
                        =np.log2(data_lbl.loc[:,('Ni%scut' % type_NorS)].astype('float'))
                        data_lbl.loc[(data_lbl.loc[:,('Ni%scutlog' % type_NorS)]==-np.inf) \
                                     | (data_lbl.loc[:,('Ni%scutlog' % type_NorS)]==np.inf),('Ni%scutlog' % type_NorS)]=np.nan
                    
                    data_lbl.loc[:,'refi']=mutids_converter(data_lbl.loc[:,'mutids'],'refi',type_form)
                    data_lbl.loc[:,'ref']=mutids_converter(data_lbl.loc[:,'mutids'],'ref',type_form)
                    data_lbl.loc[:,'mut']=mutids_converter(data_lbl.loc[:,'mutids'],'mut',type_form)

                    sbam_fh=lbl_mat_mut_cds_fh.replace('.mat_mut_cds','')
                    if not exists(prj_dh+"/data_coverage"):
                        try:
                            makedirs(prj_dh+"/data_coverage")
                        except :
                            logging.warning("race error data_coverage")                    
                    if not exists(prj_dh+"/data_lbl/"+type_form):
                        try:
                            makedirs(prj_dh+"/data_lbl/"+type_form)
                        except :
                            logging.warning("race error data_lbl")
                    if exists(sbam_fh):
                        depth_ref_fh="%s/data_coverage/%s.depth_ref" % (prj_dh,lbli)
                        depth_ref=getdepth_ref(sbam_fh,fsta_fh,cctmr=cctmr,data_out_fh=depth_ref_fh)
                        if 'refi' in depth_ref:
                            depth_ref=depth_ref.set_index('refi')
                        if 'refi' in data_lbl:
                            data_lbl=data_lbl.set_index('refi')
                        data_lbl=data_lbl.join(depth_ref)
                        data_lbl=data_lbl.reset_index()
                    #clip ends
                    if not clips is None:
                        cols=[col for col in data_lbl.columns if not (('mut' in col) or ('ref' in col))]
                        # print cols
                        data_lbl.loc[(data_lbl.loc[:,'refi']<clips[0]),cols]=np.nan
                        data_lbl.loc[(data_lbl.loc[:,'refi']>clips[1]),cols]=np.nan
                    data_lbl.reset_index().to_csv('%s/data_lbl/%s/%s' % (prj_dh,type_form,str(lbli)),index=False)
                else:
                    logging.error("len(data_lbl)/reflen is %d instead of 21 or 64" % (len(data_lbl)/(reflen/3)))
            else :
                logging.info("already processed: %s" % (str(lbli)))
        if not exists(prj_dh+"/data_mutmat"):
            try:
                makedirs(prj_dh+"/data_mutmat")
            except :
                logging.warning("race error data_mutmat")
        lbl_mat_cds_out_fh='%s/data_mutmat/%s' % (prj_dh,basename(lbl_mat_mut_cds_fh))
        lbl_mat_cds.to_csv(lbl_mat_cds_out_fh)
    else :
        logging.warning("can not find lbl_mat_mut_cds_fh : %s" % (lbl_mat_mut_cds_fh))

def data_lbl2data_fit(unsel_lbl,sel_lbl,norm_type,prj_dh,cctmr,lbls,fsta_fh):
    """
    This estimates Fitness (data_fit) from mutation data (data_lbl) in selcted and unselected samples.  
    
    :param unsel_lbl: name of input(unselected) sample.
    :param sel_lbl: name of selected sample.
    :param norm_type: default: wrt 'wild' type ['wild','syn','none'].
    :param prj_dh: path to project directory.
    :lbls: dataframe of `lbls` configuration file.
    """
    fit_lbl=sel_lbl+"_WRT_"+unsel_lbl
    for type_form in mut_types_form : # cds OR aas
        if not exists("%s/data_fit/%s/%s" % (prj_dh,type_form,fit_lbl)):             
            if  (exists("%s/data_lbl/%s/%s" % (prj_dh,type_form,sel_lbl))) \
            and (exists("%s/data_lbl/%s/%s" % (prj_dh,type_form,unsel_lbl))) : 
                logging.info("processing: data_fit/%s/%s" % (type_form,fit_lbl))        
                unsel_data=pd.read_csv(('%s/data_lbl/%s/%s' % (prj_dh,type_form,unsel_lbl)))
                sel_data  =pd.read_csv(('%s/data_lbl/%s/%s' % (prj_dh,type_form,  sel_lbl)))
                data_fit=concat_cols(unsel_data,sel_data,'mutids',
                        ['mut','ref','NiAcutlog',"NiScutlog",'depth_ref'],
                                     ['NiAcutlog',"NiScutlog",'depth_ref'],
                        'unsel','sel')
                data_fit.loc[:,'NiAunsel']=data_fit.loc[:,'NiAcutlogunsel']
                data_fit.loc[:,'NiAsel']=data_fit.loc[:,'NiAcutlogsel']
                data_fit.loc[:,'NiSunsel']=data_fit.loc[:,'NiScutlogunsel']
                data_fit.loc[:,'NiSsel']=data_fit.loc[:,'NiScutlogsel']                
                data_fit.loc[:,'FCA']=data_fit.loc[:,'NiAsel']-data_fit.loc[:,'NiAunsel'] # fold change all
                data_fit.loc[:,'FCS']=data_fit.loc[(data_fit.loc[:,'mut']==data_fit.loc[:,'ref']),'FCA']
                data_fit.loc[:,'FC_depth_ref']=\
                np.log2(data_fit.loc[:,'depth_refsel'])\
                -np.log2(data_fit.loc[:,'depth_refunsel'])
                data_fit.loc[:,'FCA_norm']=np.nan
                data_fit=data_fit.reset_index()
                if 'wild' in norm_type:                            
                    # data_fit['FiA'],
                    # data_fit.loc[(data_fit.loc[:,'mut']==data_fit.loc[:,'ref']),'FCW']\
                    # =data_fit2norm_wrt_wild(unsel_lbl,sel_lbl,type_form,data_fit['FCA'],cctmr,lbls,fsta_fh)
                    # =data_fit2norm_wrt_wild(unsel_lbl,sel_lbl,type_form,data_fit['FCA'],cctmr,lbls)
                    data_fit.loc[:,'FCA_norm']=data_fit.loc[:,'FCA']-data_fit.loc[:,'FC_depth_ref']
                elif 'syn' in norm_type:
                    try:
                        gauss_mean, gauss_sdev = fit_gauss_params(data_fit.loc[:,'FCS'])
                        data_fit.loc[:,'FCA_norm']=(data_fit.loc[:,'FCA']-gauss_mean)/gauss_sdev
                    except:
                        logging.info("norm wrt syn: excepted: data_fit/%s/%s" % (type_form,fit_lbl))
                elif norm_type == 'none':
                    data_fit.loc[:,'FCA_norm']=data_fit.loc[:,'FCA']
                data_fit.loc[:,'FiS']=data_fit.loc[(data_fit.loc[:,'mut']==data_fit.loc[:,'ref']),'FCA_norm']
                # print data_fit.head()
                data_fit=rescale_fitnessbysynonymous(data_fit,col_fit="FCA_norm",col_fit_rescaled="FiA")
                data_fit=class_fit(data_fit)
                if not exists(prj_dh+"/data_fit/"+type_form):
                    try:
                        makedirs(prj_dh+"/data_fit/"+type_form)
                    except:
                        logging.info("race error /data_fit/")
                data_fit.reset_index().to_csv('%s/data_fit/%s/%s' % (prj_dh,type_form,fit_lbl),index=False)
            else :
                logging.warning("data_lbl not present: %s/%s or %s/%s" % (type_form,sel_lbl,type_form,unsel_lbl))
        else :
            logging.info("already processed: %s" % (fit_lbl))

def class_fit(data_fit_df,zscore=False): #column of the data_fit
    """
    This classifies the fitness of mutants into beneficial, neutral or, deleterious.
    
    :param data_fit_df: dataframe of `data_fit`.
    :returns data_fit_df: classes of fitness written in 'class-fit' column based on values in column 'FiA'. 
    """
    if not zscore:
        data_fit_df.loc[data_fit_df.loc[:,'FiA']>0,'class_fit']='beneficial'
        data_fit_df.loc[data_fit_df.loc[:,'FiA']<0,'class_fit']='deleterious'
        data_fit_df.loc[data_fit_df.loc[:,'FiA']==0,'class_fit']='neutral'
    else:
        data_fit_df.loc[data_fit_df.loc[:,'FiA']>=+2,    'class_fit']="beneficial"
        data_fit_df.loc[((data_fit_df.loc[:,'FiA']>-2) & (data_fit_df.loc[:,'FiA']<+2)),'class_fit']="neutral"
        data_fit_df.loc[data_fit_df.loc[:,'FiA']<=-2,    'class_fit']="deleterious"
    return data_fit_df

def rescale_fitnessbysynonymous(data_fit,col_fit="FCA",col_fit_rescaled="FiA"):
    if col_fit_rescaled in data_fit.columns:
        col_fit_rescaled_ori=col_fit_rescaled
        col_fit_rescaled    ="tmp"
    if not "refrefi" in data_fit.columns:
        from dms2dfe.lib.io_nums import str2num
        data_fit.loc[:,'refrefi']=\
        [("%s%03d" % (mutid[0],str2num(mutid))) for mutid in data_fit.loc[:,"mutids"].tolist()]        
    for refrefi in data_fit.loc[:,"refrefi"].unique():
        subset=data_fit.loc[data_fit.loc[:,"refrefi"]==refrefi,:]
        FiW=float(subset.loc[subset.loc[:,"mut"]==subset.loc[:,"ref"],col_fit])
        for subseti in subset.index.values:
            data_fit.loc[subseti,col_fit_rescaled]=data_fit.loc[subseti,col_fit]-FiW
    if "tmp" in data_fit.columns:
        data_fit.loc[:,col_fit_rescaled_ori]=data_fit.loc[:,"tmp"]
        data_fit=data_fit.drop("tmp",axis=1)
    return data_fit

def class_comparison(data_comparison):
    """
    This classifies differences in fitness i.e. relative fitness into positive, negative or robust categories. 
    
    :param data_comparison: dataframe with `data_comparison`. 
    :returns data_comparison: dataframe with `class__comparison` added according to fitness levels in input and selected samples in `data_comparison`
    """
    if not (all(pd.isnull(data_comparison.loc[:,"class_fit_test"]))\
    or all(pd.isnull(data_comparison.loc[:,"class_fit_test"]))):
        data_comparison.loc[ ((data_comparison.loc[:,"class_fit_test"]=='beneficial')   & (data_comparison.loc[:,"class_fit_ctrl"]=='deleterious')), 'class_comparison']="positive"
        data_comparison.loc[ ((data_comparison.loc[:,"class_fit_test"]=='deleterious')  & (data_comparison.loc[:,"class_fit_ctrl"]=='beneficial')), 'class_comparison']="negative"
        data_comparison.loc[ (data_comparison.loc[:,"class_fit_test"]==data_comparison.loc[:,"class_fit_ctrl"]) , 'class_comparison']="robust"
    return data_comparison

def data_fit2data_comparison(lbl_ctrl,lbl_test,prj_dh):
    """
    This converts the Fitness values to Relative fitness among test and control fed as input.
    
    .. code-block:: text
    
        class fit:           beneficial, neutral, deleterious
        class comparison:    killed, survived
        class comparison:    positive, negative, robust
        
    :param lbl_ctrl: name of control sample.
    :param lbl_test: name of test sample.
    :param prj_dh: path to project directory.
    """    
    logging.info("processing: ctrl, test : %s %s" % (lbl_ctrl,lbl_test))
    for type_form in mut_types_form : # get aas or cds
        data_fit_ctrl_fhs=glob("%s/data_fit/%s/%s_WRT*" % (prj_dh,type_form,lbl_ctrl))
        data_fit_ctrl_keys=[basename(fh) for fh in data_fit_ctrl_fhs]
        data_fit_test_fhs=glob("%s/data_fit/%s/%s_WRT*" % (prj_dh,type_form,lbl_test))
        data_fit_test_keys=[basename(fh) for fh in data_fit_test_fhs]
        if (data_fit_ctrl_keys and data_fit_test_keys):
            for ctrli in data_fit_ctrl_keys :
                for testi in data_fit_test_keys :       
                    data_comparison_fh='%s/data_comparison/%s/%s_VERSUS_%s' % (prj_dh,type_form,testi,ctrli)
                    if data_comparison_fh.count("inferred")!=1:
                        # print data_comparison_fh
                        data_fit_ctrl=pd.read_csv("%s/data_fit/%s/%s" % (prj_dh,type_form,ctrli))
                        data_fit_test=pd.read_csv("%s/data_fit/%s/%s" % (prj_dh,type_form,testi))                    
                        data_comparison=concat_cols(data_fit_test,data_fit_ctrl,'mutids',
                                ['mut','ref','FiA',"class_fit"],['FiA',"class_fit"],
                                '_test','_ctrl')
                        
                        data_comparison.loc[:,"Fi_ctrl"]=data_comparison.loc[:,"FiA_ctrl"]
                        data_comparison.loc[:,"Fi_test"]=data_comparison.loc[:,"FiA_test"]
    #                     data_comparison.to_csv("test_comparison")
                        data_comparison=class_comparison(data_comparison) # get class fit rel
                        if not exists('%s/data_comparison/%s' % (prj_dh,type_form)):
                            try:
                                makedirs('%s/data_comparison/%s' % (prj_dh,type_form))
                            except:
                                logging.info("race error data_comparison")
                        data_comparison.reset_index().to_csv(data_comparison_fh,index=False) # store                    
        else:
            logging.warning("do not exist: data_fit/%s/%s & data_fit/%s/%s" % (type_form,lbl_ctrl,type_form,lbl_test))
