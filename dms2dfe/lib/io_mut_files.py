#!usr/bin/python

# Copyright 2016, Rohan Dandage <rraadd_8@hotmail.com,rohan@igib.in>
# This program is distributed under General Public License v. 3.  

"""
================================
``io_mut_files``
================================
"""
import sys
from os.path import splitext,exists,basename
from os import makedirs,stat
import pandas as pd
import numpy as np
from glob import glob
import pysam
import logging
logging.basicConfig(format='[%(asctime)s] %(levelname)s\tfrom %(filename)s in %(funcName)s(..): %(message)s',level=logging.DEBUG) # filename=cfg_xls_fh+'.log'
from dms2dfe.lib.fit_curve import fit_gauss_params
from dms2dfe.lib.global_vars import aas_21,cds_64,mut_types_form,mut_types_NorS
from dms2dfe.lib.convert_seq import cds2aas

def makemutids(data_lbl,refis):
    """
    This makes the mutation ids eg. M001T or A024T from mutation matrix.
    
    :param data_lbl: dataframe with columns of reference (`ref`) and mutation (`mut`).
    :param refis: index of reference amino acids/ codons.
    :returns mutids: list of mutation ids.
    """
    mutids=[]
    data_lbl=data_lbl.reset_index()
    ref_len=len(refis)
    for muti in range(len(data_lbl)/ref_len):
        refii=0
        for refi in refis:
            mutids.append("%s%03d%s" % (data_lbl.loc[muti*ref_len+refii,'ref'],refi,data_lbl.loc[muti*ref_len+refii,'mut']))
            refii+=1
    if len(data_lbl)!=len(mutids):
        logging.error("len(data_lbl)!=len(mutids) bcz %d != %d" % (len(data_lbl),len(mutids)))
        sys.exit()
    return mutids    
    
def mat_cds2mat_aas(mat_cds,host) :
    """
    This converts mutation matrix of codon to that of amino acid.
    
    :param mat_cds: codon level mutation matrix.
    :param host: name of host organism for choosing codon table. [coli | yeast | sapiens].
    :returns mat_aas: amino acid level mutation matrix.
    """
    #aas_21=["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","X"] #for indexing
    aas_wt=[cds2aas(cd,host) for cd in list(mat_cds.index)]
    aas_64=[cds2aas(cd,host) for cd in list(mat_cds.columns)]
    mat_aas_tmp=pd.DataFrame(np.float64(mat_cds.loc[:,:]),columns=aas_64,index=aas_wt) # copt cds dataframe as aas
    mat_aas    =pd.DataFrame(columns=aas_21,index=aas_wt)
    for aai in aas_21 : 
        seriesORdf= len(np.shape(mat_aas_tmp.loc[:,'%s' % aai]))
        if seriesORdf >1 :
            mat_aas.loc[:,'%s' % aai]=mat_aas_tmp.loc[:,'%s' % aai].sum(axis=1)
        else :
            mat_aas.loc[:,'%s' % aai]=mat_aas_tmp.loc[:,'%s' % aai]
    mat_aas.columns.name='mut'
    mat_aas.index.name='ref'
    # mat_aas.loc[:,'refi']=mat_cds.loc[:,'refi']
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
        
def getwildtypecov(lbl,lbls,cctmr=None):
    """
    This gets the codon level coverage of .fastq files.
    
    :param lbl: name of sample.
    :param lbls: dataframe with configuration `lbls`.
    :param cctmr: tuple with (1-based) boundaries of concatamers.
    :returns wt: dataframe with wild-type coverage.
    """
    fhs=glob(lbls.loc[lbl,'fhs_1']+"*")
    sbam_fh=[fh for fh in fhs if ((".s.bam" in fh) and (not "bam." in fh))][0]
    lbl_mat_mut_cds_fh=[fh for fh in fhs if "bam.mat_mut_cds" in fh][0]

    samfile = pysam.Samfile(sbam_fh, "rb" )

    cov=pd.DataFrame(columns=['cov'])
    for pileupcol in samfile.pileup("aph_wt_nt_cctmr", 0,1596,max_depth=10000000):
        cov.loc[pileupcol.pos, 'cov']=pileupcol.n
    samfile.close()

    if cctmr!=None:
        cov_cctmr1=cov.iloc[(cctmr[0][0]-1)*3:(cctmr[0][1]-1)*3,:]
        cov_cctmr2=cov.iloc[(cctmr[1][0]-1)*3:(cctmr[1][1]-1)*3,:]
        cov=cov_cctmr1.reset_index(drop=True).astype(int)+cov_cctmr2.reset_index(drop=True).astype(int)

    cov_per_cds=cov.groupby(np.array([cov.index.values, cov.index.values,cov.index.values]).T.flatten()[:len(cov.index.values)]).mean()

    # minus the mutants
    lbl_mat_mut_cds=pd.read_csv(lbl_mat_mut_cds_fh)
    if cctmr!=None:
        lbl_mat_mut_cds_cctmr1=lbl_mat_mut_cds.iloc[(cctmr[0][0]-1):(cctmr[0][1]-1),:]
        lbl_mat_mut_cds_cctmr2=lbl_mat_mut_cds.iloc[(cctmr[1][0]-1):(cctmr[1][1]-1),:]
        lbl_mat_mut_cds=lbl_mat_mut_cds_cctmr1.fillna(0).set_index("ref_cd",drop=True)+lbl_mat_mut_cds_cctmr2.fillna(0).set_index("ref_cd",drop=True)
    del lbl_mat_mut_cds['Unnamed: 0']

    cov_mut=pd.DataFrame({"cov":lbl_mat_mut_cds.T.sum().reset_index(drop=True)})

    wt=cov_per_cds - cov_mut 
    return wt

def data_fit2norm_wrt_wild(unsel_lbl,sel_lbl,type_form,FCA,cctmr,lbls):
    """
    This normalises Fold changes (FC) to normalised fitness values (Fi)
    
    Gets fhs of sbam files and using pysam, gets the coverage per nt
    bins it into 3 nt bins i.e. codons.
    substract mat mut cd sum,
    
    get log2 ratio of vectors for both sel to unsel (FCW)
    
    .. math::

        DMut = (Mut selected /WT selected )/(Mut input /WT input )
        Fi = FC - FCW
    
    :param unsel_lbl: name of input sample.
    :param sel_lbl: name of selected sample.
    :param type_form: type of mutations codon (`cds`) or amino acid (`aas`) form.
    :param FCA: Fold change(FC) values of all(A) mutations.
    :param cctmr: tuple with (1-based) boundaries of concatamers.
    :param lbls: dataframe with configuration `lbls`.
    :returns FiA:  normalised fitness values (Fi)
    """
    FCW=pd.DataFrame()
    FCW.loc[:,'FCW']=np.log(getwildtypecov(sel_lbl,lbls,cctmr)/getwildtypecov(unsel_lbl,lbls,cctmr))

    FiA=pd.DataFrame()
    if type_form=='aas':
        iterations=21
    if type_form=='cds':
        iterations=64
        
    if len(FCA)/len(FCW)==iterations:
        for i in range(iterations):
#             print "%d %d %d" % (i*len(FCW),i*len(FCW)+len(FCW)-1, i*len(FCW)+len(FCW)-i*len(FCW))
            FCA_frag=FCA.iloc[i*len(FCW):i*len(FCW)+len(FCW)]
            FCA_frag=FCA_frag.reset_index()
            FiA_frag=FCA_frag["FCA"]-FCW["FCW"]
            FiA=pd.concat([FiA,FiA_frag],ignore_index=True)
    else: 
        logging.error("len(FCA)/len(FCW) is %d instead of %d" % (len(FCA)/len(FCW),iterations))
    if len(FiA)!=len(FCA):
        logging.error("len(FiA) is %d instead of %d" % (len(FiA),len(FCA)))
    return FiA,FCW    
    
def class_fit(data_fit_df): #column of the data_fit
    """
    This classifies the fitness of mutants into beneficial, neutral or, deleterious.
    
    :param data_fit_df: dataframe of `data_fit`.
    :returns data_fit_df: classes of fitness written in 'class-fit' column based on values in column 'FiA'. 
    """
    data_fit_df.loc[data_fit_df.loc[:,'FiA']>=+2,    'class_fit']="beneficial"
    data_fit_df.loc[((data_fit_df.loc[:,'FiA']>-2) & (data_fit_df.loc[:,'FiA']<+2)),'class_fit']="neutral"
    data_fit_df.loc[data_fit_df.loc[:,'FiA']<=-2,    'class_fit']="deleterious"
    return data_fit_df
                                    

                
def getusable_lbls_list(prj_dh):
    """
    This detects the samples that can be processed.
    
    :param prj_dh: path to project directory.
    :returns lbls_list: list of names of samples that can be processed.
    """
    lbls=pd.read_csv(prj_dh+'/cfg/lbls')
    lbls=lbls.set_index('varname')
    lbls_list=[]
    #data_lbl cols: NiA mutids NiS NiN NiNcut NiNcutlog NiScut NiScutlog NiAcut NiAcutlog    
    for lbli,lbl in lbls.iterrows() :
        if (not exists("%s/data_lbl/%s/%s" % (prj_dh,'aas',str(lbli)))) \
        and (not exists("%s/data_lbl/%s/%s" % (prj_dh,'cds',str(lbli)))):
            fh_1=str(lbl['fhs_1'])
            # print fh_1
            lbl_mat_mut_cds_fh=[fh for fh in glob(fh_1+"*") if '.mat_mut_cds' in fh]
            if len(lbl_mat_mut_cds_fh)!=0:
                lbl_mat_mut_cds_fh=lbl_mat_mut_cds_fh[0]
                lbls_list.append([lbli,lbl_mat_mut_cds_fh])
            else :
                logging.warning("can not find: %s" % basename(fh_1))
        else:
            logging.info("already processed: %s" % (str(lbli)))
    return lbls_list

def getusable_fits_list(prj_dh):
    """
    This gets the list of samples that can be processed for fitness estimations.
    
    :param prj_dh: path to project directory.
    :returns fits_pairs_list: list of tuples with names of input and selected samples.
    """
    if not exists('%s/cfg/fit'% (prj_dh)):
        logging.warning("ana3_mutmat2fit : getusable_fits_list : not fits in cfg/fit")
        return []
    else:
        fits=pd.read_csv(prj_dh+'/cfg/fit')
        if "Unnamed: 0" in fits.columns:
            fits=fits.drop("Unnamed: 0", axis=1)
        fits=fits.set_index('unsel')
        fits_pairs_list=[]
        for unsel_lbl,sels in fits.iterrows() :
            sels=sels[~sels.isnull()]
            for sel_lbl in sels :
                fit_lbl=sel_lbl+"_WRT_"+unsel_lbl
                if (not exists("%s/data_fit/%s/%s" % (prj_dh,'aas',fit_lbl))) \
                and (not exists("%s/data_fit/%s/%s" % (prj_dh,'cds',fit_lbl))):             
                    if  ((exists("%s/data_lbl/%s/%s" % (prj_dh,'aas',sel_lbl))) \
                    and (exists("%s/data_lbl/%s/%s" % (prj_dh,'aas',unsel_lbl)))) or \
                        ((exists("%s/data_lbl/%s/%s" % (prj_dh,'cds',sel_lbl))) \
                    and (exists("%s/data_lbl/%s/%s" % (prj_dh,'cds',unsel_lbl)))) : 
                        fits_pairs_list.append([unsel_lbl,sel_lbl])
                    else :
                        logging.warning("data_lbl not present : %s or %s" % (sel_lbl,unsel_lbl))
                else :
                    logging.info("already processed: %s" % (fit_lbl))
        return fits_pairs_list


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
                    else:
                        logging.info("already processed %s" % avg_lbl)
    else:
        logging.warning("skipping repli2data_lbl_avg")

def mut_mat_cds2data_lbl(lbli,lbl_mat_mut_cds_fh,host,prj_dh,fsta_seqlen,cctmr,Ni_cutoff):    
    """
    This function converts mutation matrix (produced from .mat) file into data_lbl format.
    
    :param lbli: name of the sample.
    :param lbl_mat_mut_cds_fh: path to codon level mutation matrix (.mat_mut_cds).
    :param host : name of host organism for choosing codon table. (coli | yeast | sapiens).
    :param prj_dh: path to project directory.
    :param fsta_seqlen: nucleotide level length of the reference sequence.
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
            fsta_seqlen=(cctmr[0][1]-1)*3
            lbl_mat_cds=collate_cctmr(lbl_mat_cds,cctmr)
        lbl_mat_aas=mat_cds2mat_aas(lbl_mat_cds,host)# convert to aa data
        for type_form in mut_types_form : # get aas or cds  
            if not exists("%s/data_lbl/%s/%s" % (prj_dh,type_form,str(lbli))):
                if type_form=="aas":
                    data_lbl=pd.DataFrame(lbl_mat_aas.unstack())
                elif type_form=="cds":
                    data_lbl=pd.DataFrame(lbl_mat_cds.drop("refi",axis=1).unstack())
                if (len(data_lbl)/(fsta_seqlen/3)==21) or (len(data_lbl)/(fsta_seqlen/3)==64):
                    data_lbl.columns=["NiA"]  # rename col Ni  
                    data_lbl.loc[:,'mutids']=makemutids(data_lbl,lbl_mat_cds.loc[:,"refi"]) # add mutids
                    data_lbl=pd.concat([data_lbl,getNS(data_lbl)], axis=1)
                    for type_NorS in mut_types_NorS : # N OR S  
                        data_lbl.loc[:,('Ni%scut' % type_NorS)]=data_lbl.loc[:,'Ni%s' % type_NorS]
                        data_lbl.loc[data_lbl.loc[:,('Ni%scut' % type_NorS)]<Ni_cutoff,('Ni%scut' % type_NorS)]=np.nan \
                        #Ni_cutoff=8 varscan default
                        data_lbl.loc[:,('Ni%scutlog' % type_NorS)]=np.log2(data_lbl.loc[:,('Ni%scut' % type_NorS)].astype('float'))
                        data_lbl.loc[(data_lbl.loc[:,('Ni%scutlog' % type_NorS)]==-np.inf) \
                                     | (data_lbl.loc[:,('Ni%scutlog' % type_NorS)]==np.inf),('Ni%scutlog' % type_NorS)]=np.nan
                    if not exists(prj_dh+"/data_lbl/"+type_form):
                        try:
                            makedirs(prj_dh+"/data_lbl/"+type_form)
                        except :
                            logging.warning("race error data_lbl")
                    data_lbl.reset_index().to_csv('%s/data_lbl/%s/%s' % (prj_dh,type_form,str(lbli)),index=False)
                else:
                    logging.error("len(data_lbl)/ref_len is %d instead of 21 or 64" % (len(data_lbl)/(fsta_seqlen/3)))
            else :
                logging.info("already processed: %s" % (str(lbli)))
    else :
        logging.warning("can not find lbl_mat_mut_cds_fh : %s" % (lbl_mat_mut_cds_fh))

def data_lbl2data_fit(unsel_lbl,sel_lbl,norm_type,prj_dh,cctmr,lbls):
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
                data_fit=pd.DataFrame(columns=['mutids','NiAunsel','NiAsel','NiSunsel','NiSsel','FCA','FCS','FiA','FiS'], index=unsel_data.index)
                data_fit.loc[:,'mut']=unsel_data.loc[:,'mut']
                data_fit.loc[:,'ref']=unsel_data.loc[:,'ref']
                data_fit.loc[:,'mutids']=unsel_data.loc[:,'mutids']
                data_fit.loc[:,'NiAunsel']=unsel_data.loc[:,'NiAcutlog']
                data_fit.loc[:,'NiAsel']=sel_data.loc[:,'NiAcutlog']
                data_fit.loc[:,'NiSunsel']=unsel_data.loc[:,'NiScutlog']
                data_fit.loc[:,'NiSsel']=sel_data.loc[:,'NiScutlog']                
                data_fit.loc[:,'FCA']=data_fit.loc[:,'NiAsel']-data_fit.loc[:,'NiAunsel'] # fold change all
                data_fit.loc[:,'FCS']=data_fit.loc[(data_fit.loc[:,'mut']==data_fit.loc[:,'ref']),'FCA']
                if 'wild' in norm_type:                            data_fit['FiA'],data_fit.loc[(data_fit.loc[:,'mut']==data_fit.loc[:,'ref']),'FCW']=data_fit2norm_wrt_wild(unsel_lbl, \
                                                                                                    sel_lbl,type_form, \
                                                                                                    data_fit['FCA'], \
                                                                                                    cctmr,lbls)
                elif 'syn' in norm_type:
                    try:
                        gauss_mean, gauss_sdev = fit_gauss_params(data_fit.loc[:,'FCS'])
                        data_fit.loc[:,'FiA']=(data_fit.loc[:,'FCA']-gauss_mean)/gauss_sdev
                    except:
                        logging.info("norm wrt syn: excepted: data_fit/%s/%s" % (type_form,fit_lbl))
                elif norm_type == 'none':
                    data_fit.loc[:,'FiA']=data_fit.loc[:,'FCA']
                data_fit.loc[:,'FiS']=data_fit.loc[(data_fit.loc[:,'mut']==data_fit.loc[:,'ref']),'FiA']
                data_fit=class_fit(data_fit)
                if not exists(prj_dh+"/data_fit/"+type_form):
                    try:
                        makedirs(prj_dh+"/data_fit/"+type_form)
                    except:
                        logging.info("race error /data_fit/")
                data_fit.reset_index().to_csv('%s/data_fit/%s/%s' % (prj_dh,type_form,fit_lbl),index=False)
            else :
                logging.warning("data_lbl not present: %s or %s" % (sel_lbl,unsel_lbl))
        else :
            logging.info("already processed: %s" % (fit_lbl))

def rescale_fitnessbysynonymous(data_fit,col_fit="FiA",col_fit_rescaled="FiA rescaled"):
    for refrefi in data_fit.loc[:,"refrefi"].unique():
        subset=data_fit.loc[data_fit.loc[:,"refrefi"]==refrefi,:]
        FiW=float(subset.loc[subset.loc[:,"mut"]==subset.loc[:,"ref"],col_fit])
        for subseti in subset.index.values:
            data_fit.loc[subseti,col_fit_rescaled]=data_fit.loc[subseti,col_fit]-FiW
    return data_fit          

def class_comparison(data_comparison):
    """
    This classifies differences in fitness i.e. relative fitness into positive, negative or robust categories. 
    
    :param data_comparison: dataframe with `data_comparison`. 
    :returns data_comparison: dataframe with `class__comparison` added according to fitness levels in input and selected samples in `data_comparison`
    """
    if not (all(pd.isnull(data_comparison.loc[:,"class_fit_test"]))\
    or all(pd.isnull(data_comparison.loc[:,"class_fit_test"]))):
        data_comparison.loc[ ((data_comparison.loc[:,"class_fit_test"]=='beneficial')   & ~(data_comparison.loc[:,"class_fit_ctrl"]=='beneficial')) | \
                             ((data_comparison.loc[:,"class_fit_test"]=='neutral')      & (data_comparison.loc[:,"class_fit_ctrl"]=='deleterious')), 'class_comparison']="positive"
        data_comparison.loc[ ((data_comparison.loc[:,"class_fit_test"]=='deleterious')  & ~(data_comparison.loc[:,"class_fit_ctrl"]=='deleterious')) | \
                             ((data_comparison.loc[:,"class_fit_test"]=='neutral')      & (data_comparison.loc[:,"class_fit_ctrl"]=='beneficial')), 'class_comparison']="negative"
        data_comparison.loc[ (data_comparison.loc[:,"class_fit_test"]==data_comparison.loc[:,"class_fit_ctrl"]) , 'class_comparison']="robust"

        data_comparison.loc[(                  data_comparison.loc[:,"class_fit_ctrl"].isnull() & \
                            np.logical_not(data_comparison.loc[:,"class_fit_test"].isnull())), 'class_comparison']="survived"
        data_comparison.loc[(np.logical_not(data_comparison.loc[:,"class_fit_ctrl"].isnull()) & \
                                            data_comparison.loc[:,"class_fit_test"].isnull()), 'class_comparison']="killed"
        data_comparison.loc[(data_comparison.loc[:,"class_fit_ctrl"].isnull() & data_comparison.loc[:,"class_fit_test"].isnull()), 'class_comparison']=np.nan
    return data_comparison

def getusable_comparison_list(prj_dh):
    """
    This converts the table of tests and controls in configuration file into tuples of test and control.
    
    :param prj_dh: path to project directory.
    """
    comparisons=pd.read_csv(prj_dh+'/cfg/comparison')
    comparisons=comparisons.set_index('ctrl')
    comparison_list=[]
    for ctrl,row in comparisons.iterrows() :
        row=row[~row.isnull()]
        for test in row[0:] :
            comparison_list.append([ctrl,test])
    return  comparison_list  

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
                    data_fit_ctrl=pd.read_csv("%s/data_fit/%s/%s" % (prj_dh,type_form,ctrli))
                    data_fit_test=pd.read_csv("%s/data_fit/%s/%s" % (prj_dh,type_form,testi))                    
#                     data_fit_ctrl.to_csv("test_ctrl")
#                     data_fit_test.to_csv("test_test")
#                     data_comparison=pd.DataFrame({'mutids':data_fit_ctrl.loc[:,"mutids"], \
#                                                 'Fi_ctrl':data_fit_ctrl.loc[:,"FiA"],\
#                                                 'class_fit_ctrl':data_fit_ctrl.loc[:,"class_fit"], \
#                                                 'Fi_test':data_fit_test.loc[:,"FiA"],\
#                                                 'class_fit_test':data_fit_test.loc[:,"class_fit"]})
                    data_comparison=pd.DataFrame()
                    data_comparison.loc[:,"mutids"]=data_fit_ctrl.loc[:,"mutids"]
                    data_comparison.loc[:,"mut"]=data_fit_ctrl.loc[:,"mut"]
                    data_comparison.loc[:,"ref"]=data_fit_ctrl.loc[:,"ref"]
                    data_comparison.loc[:,"Fi_ctrl"]=data_fit_ctrl.loc[:,"FiA"]
                    data_comparison.loc[:,"class_fit_ctrl"]=data_fit_ctrl.loc[:,"class_fit"]
                    data_comparison.loc[:,"Fi_test"]=data_fit_test.loc[:,"FiA"]
                    data_comparison.loc[:,"class_fit_test"]=data_fit_test.loc[:,"class_fit"]
#                     data_comparison.to_csv("test_comparison")
                    data_comparison=class_comparison(data_comparison) # get class fit rel
                    data_comparison_fh='%s/data_comparison/%s/%s_VERSUS_%s' % (prj_dh,type_form,testi,ctrli)
                    if not exists('%s/data_comparison/%s' % (prj_dh,type_form)):
                        try:
                            makedirs('%s/data_comparison/%s' % (prj_dh,type_form))
                        except:
                            logging.info("race error data_comparison")
                    data_comparison.reset_index().to_csv(data_comparison_fh,index=False) # store                    
        else:
            logging.warning("do not exist: data_fit/%s/%s & data_fit/%s/%s" % (type_form,lbl_ctrl,type_form,lbl_test))