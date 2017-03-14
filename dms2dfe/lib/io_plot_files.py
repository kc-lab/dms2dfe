#!usr/bin/python

# Copyright 2016, Rohan Dandage <rraadd_8@hotmail.com,rohan@igib.in>
# This program is distributed under General Public License v. 3.  

"""
================================
``io_plot_files``
================================
"""
import sys
from os import makedirs,stat
from os.path import splitext, join, exists, isdir,basename,abspath,dirname
import pandas as pd
import numpy as np
import subprocess
from glob import glob
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from Bio import SeqIO
from pychimera.pychimera import guess_chimera_path
import logging
logging.basicConfig(format='[%(asctime)s] %(levelname)s\tfrom %(filename)s in %(funcName)s(..): %(message)s',level=logging.DEBUG) # filename=cfg_xls_fh+'.log'
from dms2dfe import configure
from dms2dfe.lib.plot_mut_data import plot_data_lbl_repli,plot_data_fit_scatter,plot_data_fit_dfe,plot_data_fit_heatmap,data2mut_matrix,plot_data_comparison_bar,plot_data_fit_clustermap,plot_sub_matrix,plot_cov,plot_data_comparison_violin
# from dms2dfe.lib.io_seq_files import getwildtypecov
from dms2dfe.lib.plot_pdb import vector2bfactor
from dms2dfe.lib.global_vars import mut_types_form
from dms2dfe.lib.io_dfs import denanrows
# from dms2dfe.lib.io_plot_files import plot_coverage,plot_heatmap,plot_clustermap,plot_multisca,plot_violin,plot_pies,plot_pdb

def plot_coverage(info,plot_type="coverage",type_form='aas'):    
    if not exists(info.prj_dh+"/data_coverage/aas"):
        makedirs(info.prj_dh+"/data_coverage/aas")
    lbls=pd.read_csv("%s/cfg/lbls" % info.prj_dh).set_index("varname",drop=True)
    for lbl in lbls.index.values:
        plot_fh="%s/plots/%s/%s.%s.pdf" % (info.prj_dh,type_form,lbl,plot_type)
        data_fh="%s/data_coverage/%s.depth_ref" % (info.prj_dh,lbl)        
        if not exists(plot_fh):
            if not pd.isnull(lbls.loc[lbl,'fhs_1']):
                fhs=glob(lbls.loc[lbl,'fhs_1']+"*")
                # print fhs
                if len(fhs)!=0:
                    sbam_fhs=[fh for fh in fhs if (fh.endswith(".s.bam"))]
                    if len(sbam_fhs)!=0:
                        sbam_fh=sbam_fhs[0]
                        lbl_mat_mut_cds_fh=[fh for fh in fhs if "bam.mat_mut_cds" in fh][0]
                        # print lbl_mat_mut_cds_fh
                        if exists(data_fh):
                            data_cov=pd.read_csv(data_fh)
                        data_lbl=pd.read_csv("%s/data_lbl/aas/%s" % (info.prj_dh,lbl))
                        plot_cov(data_cov,data_lbl,plot_fh=plot_fh)
                else:
                    logging.warning("can not find sequencing data to get coverage")
            else:
                logging.warning("can not find sequencing data to get coverage")

def get_data_fhs(info):
    data_fit_fhs=glob("%s/data_fit/aas/*" % info.prj_dh)
    # +glob("%s/data_fit/cds/*" % info.prj_dh)
    data_fit_fhs= [fh for fh in data_fit_fhs if "_WRT_" in basename(fh)]
    data_fit_fhs=np.unique(data_fit_fhs)
    return data_fit_fhs

def get_fhs(dh,include,exclude=None):
    fhs=glob("%s/*" % dh)
    fhs= [fh for fh in fhs if include in basename(fh)]
    if not exclude is None:
        fhs= [fh for fh in fhs if not exclude in basename(fh)]
    fhs=np.unique(fhs)
    return fhs
    
def plot_heatmap(info,data_fit_fhs=None,plot_type="heatmap"):
    if data_fit_fhs is None:
        data_fit_fhs=get_data_fhs(info)
    type_form='aas'
    data_feats=pd.read_csv(info.prj_dh+"/data_feats/aas/feats_all")

    for data_fit_fh in data_fit_fhs:
        data_fit_fn=basename(data_fit_fh)
        data_fit=pd.read_csv(data_fit_fh)
        plot_fh="%s/plots/%s/%s.%s.pdf" % (info.prj_dh,type_form,data_fit_fn,plot_type) 
        if not exists(plot_fh):
            if "_inferred" in basename(plot_fh):
                cbar_label="Fitness scores ($F_{i}$)"
            else:
                cbar_label="Log fold change ($FC_{i}$)"
            plot_data_fit_heatmap(data_fit,type_form,col='FiA',data_feats=data_feats,\
                                  plot_fh=plot_fh,cbar_label=cbar_label)
    
def plot_clustermap(info,data_fit_fhs=None,plot_type="clustermap"):
    if data_fit_fhs is None:
        data_fit_fhs=get_data_fhs(info)
    type_form='aas'
    for data_fit_fh in data_fit_fhs:
        data_fit_fn=basename(data_fit_fh)
        data_fit=pd.read_csv(data_fit_fh)
        plot_fh="%s/plots/%s/%s.%s.pdf" % (info.prj_dh,type_form,data_fit_fn,plot_type) 
        if not exists(plot_fh):
            ax=plot_data_fit_clustermap(data_fit,type_form,col='FiA',col_cluster=True)
            ax.savefig(plot_fh+".pdf",format='pdf')                        
            ax.savefig(plot_fh);plt.clf();plt.close()

from dms2dfe.lib.plot_mut_data_scatter import plot_scatter_mutilayered

def plot_multisca(info,data_fit_fhs=None,plot_type='multisca'):
    if data_fit_fhs is None:
        data_fit_fhs=get_data_fhs(info)
    data_fit_fhs=[fh for fh in data_fit_fhs if not '_inferred' in fh]    
    type_form='aas'
    for data_fit_fh in data_fit_fhs:
        data_fit_fn=basename(data_fit_fh)
        data_fit=pd.read_csv(data_fit_fh).set_index('mutids')
        plot_fh='%s/plots/aas/%s.%s.pdf' % (info.prj_dh,data_fit_fn,plot_type)
        if not exists(plot_fh):
            cols=['NiA_tran.ref','NiA_tran.sel']
            cols_labels=['$N_{i,ref}$','$N_{i,sel}$']
            data=data_fit.loc[denanrows(data_fit.loc[:,cols]).index.tolist(),:]
            data.loc[:,cols_labels[0]]=data.loc[:,cols[0]]
            data.loc[:,cols_labels[1]]=data.loc[:,cols[1]]
            m,s,p=plot_scatter_mutilayered(data,
                                           cols_labels[0],cols_labels[1],
            #                          mutids_heads=mutids_heads,
            #                          mutids_tails=mutids_tails,
                                     color_heads='b',color_tails='b',
                                     col_z_mutations='padj',
                                     zcol_threshold=0.05,
                                     repel=0.035,
                                     figsize=[6.5,2.5],#[6.375,4.5],
                                    color_dots='tails',
                                     diagonal=False,
                                        errorbars=True,
                                     plot_fh=plot_fh,
                                     space=0.1,
                                    )    
    
def plot_pdb(info,data_fit_fhs=None,plot_type="pdb"):
    if data_fit_fhs is None:
        data_fit_fhs=get_data_fhs(info)
    type_form='aas'
    for data_fit_fh in data_fit_fhs:
        data_fit_fn=basename(data_fit_fh)
        pdb_clrd_fh="%s/plots/%s/%s.pdb" % (info.prj_dh,type_form,data_fit_fn) 
        if not exists(pdb_clrd_fh):
            data_fit=pd.read_csv(data_fit_fh)
            mut_matrix=data2mut_matrix(data_fit,'FiA','mut',type_form)
            data_fit_avg=mut_matrix.mean()
            vector2bfactor(data_fit_avg,info.pdb_fh,pdb_clrd_fh)
            plot_pdb_chimera_fhs_f.write(abspath(pdb_clrd_fh)+"\n")
        else:
            logging.info("already processed: %s" % basename(pdb_clrd_fh))
    else:
        logging.info("already processed")
    plot_pdb_chimera_fhs_f.close()
    try:
        chimera_dh=guess_chimera_path()[0]
        if chimera_dh:
            std=subprocess.Popen("which glxinfo",shell=True,stdout=subprocess.PIPE)
            if std.stdout.read():
                if not stat(plot_pdb_chimera_fhs_fh).st_size == 0:
                    subprocess.call("%s/bin/chimera --silent %s/lib/plot_pdb_chimera.py" % (chimera_dh,abspath(dirname(__file__))),shell=True)
                else:
                    logging.info("already processed")  
            else:
                logging.error("skipping: pdb vizs: graphics drivers not present/configured.") 
                logging.info("To configure graphics drivers for UCSF-Chimera please install mesa-utils: sudo apt-get install mesa-utils;sudo apt-get update ")  
        else:
            logging.info("skipping pdb vizs : please install UCSF-Chimera ")      
    except:
        logging.info("skipping pdb vizs : please install UCSF-Chimera")      
    logging.shutdown()

from dms2dfe.lib.plot_mut_data_dists import plot_data_comparison_multiviolin

def plot_violin(info,data_comparison_fhs=None,plot_type='violin'):
    if data_comparison_fhs is None:
        data_comparison_fhs=get_fhs('%s/data_comparison/aas/' % info.prj_dh,
                                    include='_VERSUS_',exclude='_inferred')
    aasORcds="aas"
    col="FiA"
    ylabel="$F_{i}$"
    ylims=[-4,6]
    for fh in data_comparison_fhs:
        fn=basename(fh)
        plot_fh='%s/plots/aas/%s.%s.pdf' % (info.prj_dh,fn,plot_type)
        data_fit_fns=fn.split('_VERSUS_')[::-1]
        data_fiti_ctrl=0
        data_fit_labels=data_fit_fns
        # print plot_fh
        # print fn
        # print data_fit_fns
        plot_data_comparison_multiviolin(info.prj_dh,
                                data_fit_fns,col,
                                data_fiti_ctrl=data_fiti_ctrl,
                                aasORcds=aasORcds,
                                data_fits_labels=data_fit_labels,
                                color_test='yellow',
                                color_ctrl='lightgray',
                                figsize=[2.65,2.5],
                                color_xticks=(0,0.2,1),
                                ns=False,
                                numeric=True,
                                # ylims=ylims,
                                force=True,
                                ylabel=ylabel,
                                plot_fh=plot_fh)

from dms2dfe.lib.plot_mut_data_dists import plot_pie

def plot_pies(info,data_comparison_fhs=None,plot_type='pies_selection'):    
    if data_comparison_fhs is None:
        data_comparison_fhs=get_fhs('%s/data_comparison/aas/' % info.prj_dh,
                                    include='_VERSUS_',exclude='_inferred')
    for fh in data_comparison_fhs:
        fn=basename(fh)
        plot_fh='%s/%s.%s.pdf' % (info.prj_dh,fn,plot_type)
        data_fit=pd.read_csv(fh)
        col_classes='class_comparison'
        colors=['cyan', 'hotpink','lightgray']
        explodei='mid'
        fontsize=12.5
        flag=''
        plot_pie(data_fit,col_classes,"mutids",
                 explodei=explodei,
                 figsize=[1,1],
                 colors=colors,
                 plot_fh=plot_fh,
                 label=True,
                 fontsize=fontsize,
                )