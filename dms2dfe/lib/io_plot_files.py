#!usr/bin/python

# Copyright 2016, Rohan Dandage <rraadd_8@hotmail.com,rohan@igib.in>
# This program is distributed under General Public License v. 3.  

"""
================================
``io_plot_files``
================================
"""
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

def plot_coverage(info,plot_type="coverage"):    
    if not exists(prj_dh+"/data_coverage/aas"):
        makedirs(prj_dh+"/data_coverage/aas")
    lbls=pd.read_csv("%s/cfg/lbls" % info.prj_dh).set_index("varname",drop=True)
    type_form="aas"
    for lbl in lbls.index.values:
        plot_fh="%s/plots/%s/plot_%s_%s.pdf" % (prj_dh,type_form,plot_type,lbl)
        data_fh="%s/data_coverage/%s/%s" % (prj_dh,type_form,lbl)        
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
                        if not exists(data_fh):
                            data_cov=getwildtypecov(sbam_fh,lbl_mat_mut_cds_fh,fsta_id,fsta_seqlen,cctmr)
                            data_cov.to_csv(data_fh,index=False)
                        else: 
                            data_cov=pd.read_csv(data_fh)
                        data_lbl=pd.read_csv("%s/data_lbl/aas/%s" % (prj_dh,lbl))
                        plot_cov(data_cov,data_lbl,plot_fh=plot_fh)
                else:
                    logging.warning("can not find sequencing data to get coverage")
            else:
                logging.warning("can not find sequencing data to get coverage")

    data_fit_fhs=glob("%s/data_fit/aas/*" % prj_dh)+glob("%s/data_fit/cds/*" % prj_dh)
    data_fit_fhs= [fh for fh in data_fit_fhs if "_WRT_" in basename(fh)]
    data_fit_fhs=np.unique(data_fit_fhs)
    
def plot_heatmap(info,data_fit_fhs,plot_type="heatmap"):
    for data_fit_fh in data_fit_fhs:
        data_fit=pd.read_csv(data_fit_fh)
        plot_fh="%s/plots/%s/fig_%s_%s.pdf" % (prj_dh,type_form,plot_type,data_fiti) 
        if not exists(plot_fh):
            if "_inferred" in basename(plot_fh):
                cbar_label="Fitness scores ($F_{i}$)"
            else:
                cbar_label="Log fold change ($FC_{i}$)"
            plot_data_fit_heatmap(data_fit,type_form,col='FiA',data_feats=data_feats,\
                                  plot_fh=plot_fh,cbar_label=cbar_label)
    
def plot_clustermap(info,data_fit_fhs,plot_type="clustermap"):
    for data_fit_fh in data_fit_fhs:
        data_fit=pd.read_csv(data_fit_fh)
        plot_fh="%s/plots/%s/fig_%s_%s.pdf" % (prj_dh,type_form,plot_type,data_fiti) 
        if not exists(plot_fh):
            ax=plot_data_fit_clustermap(data_fit,type_form,col='FiA',col_cluster=True)
            ax.savefig(plot_fh+".pdf",format='pdf')                        
            ax.savefig(plot_fh);plt.clf();plt.close()

            
def plot_multisca(info,data_fit_fhs,plot_type='multisca'):
    for data_fit_fh in data_fit_fhs:
        data_fit=pd.read_csv(data_fit_fh).set_index('mutids')
        plot_fh='%s/plot_%s.pdf' % (plot_dh,plot_type)
        if not exists(plot_fh):
            print plot_fh
            cols=['NiA_tran.ref','NiA_tran.sel']
            cols_labels=['$N_{i,ref}$','$N_{i,sel}$']
            data=data.loc[denanrows(data.loc[:,cols]).index.tolist(),:]
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
    
            
def plot_pdb(info,data_fit_fhs,plot_type="pdb"):
    for data_fit_fh in data_fit_fhs:
        pdb_clrd_fh="%s/plots/%s/fig_%s_%s.pdb" % (prj_dh,type_form,plot_type,data_fiti) 
        if not exists(pdb_clrd_fh):
            data_fit=pd.read_csv(data_fit_fh)
            mut_matrix=data2mut_matrix(data_fit,'FiA','mut',type_form)
            data_fit_avg=mut_matrix.mean()
            vector2bfactor(data_fit_avg,pdb_fh,pdb_clrd_fh)
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

def plot_violin(plot_type='violin'):
        
    plot_fh='%s/plot_%s.pdf' % (plot_dh,plot_type)
    # %run ../progs/dms2dfe/dms2dfe/lib/plot_mut_data_dists.py
    prj_dh='/home/kclabws1/Documents/propro/writ/prjs/1_dms_software/data/datasets/APH2_Melnikov_et_al_2014/b170104_full_run//rlog_GLM'

    data_fit_fns=['Kan14_WRT_Bkg','Kan11_WRT_Bkg']
    data_fiti_ctrl=0
    aasORcds="aas"
    col="FiA"
    ylabel="$F_{i}$"
    ylims=[-4,6]

    data_fit_labels=[data_fit_fns2labels[s] if s in data_fit_fns2labels.keys() else s for s in data_fit_fns]

    tmp,tmp2=plot_data_comparison_multiviolin(prj_dh,data_fit_fns,col,
                                              data_fiti_ctrl=data_fiti_ctrl,
                                              aasORcds=aasORcds,
                                            data_fits_labels=data_fit_labels,#[l.replace(' + ','\n+\n') for l in data_fit_labels],ylabel=ylabel,
                                            color_test='yellow',
                                            color_ctrl='lightgray',
                                            figsize=[2.65,2.5],
                                            color_xticks=(0,0.2,1),
    #                                         plot_fh=plot_fh,
                                            ns=False,
                                            numeric=True,
                                            ylims=ylims,
                                            force=True,
                                              ylabel=ylabel,
                                              plot_fh=plot_fh
                                             )

def plot_pies():
    plot_type='plot_pies_selection'
    plot_fh='%s/plot_%s.pdf' % (plot_dh,plot_type)
    print plot_fh
    # %run ../progs/dms2dfe/dms2dfe/lib/plot_mut_data_dists.py

    col_classes='class_comparison'
    colors=['cyan', 'hotpink','lightgray']
    explodei='mid'
    fontsize=12.5
    data_fh='/home/kclabws1/Documents/propro/writ/prjs/1_dms_software/data/datasets/APH2_Melnikov_et_al_2014/b170104_full_run//rlog_GLM/data_comparison/aas/Kan11_WRT_Bkg_VERSUS_Kan14_WRT_Bkg'
    data_fit=pd.read_csv(data_fh)
    flag=''
    plot_pie(data_fit,col_classes,"mutids",
             explodei=explodei,
             figsize=[1,1],
             colors=colors,
             plot_fh=plot_fh,
             label=True,
             fontsize=fontsize,
            )