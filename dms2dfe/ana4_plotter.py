#!usr/bin/python

# Copyright 2016, Rohan Dandage <rraadd_8@hotmail.com,rohan@igib.in>
# This program is distributed under General Public License v. 3.  

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
from pychimera.pychimera import guess_chimera_path
import logging
logging.basicConfig(format='[%(asctime)s] %(levelname)s\tfrom %(filename)s in %(funcName)s(..): %(message)s',level=logging.DEBUG) # filename=cfg_xls_fh+'.log'
from dms2dfe import configure
from dms2dfe.lib.plot_mut_data import plot_data_lbl_repli,plot_data_fit_scatter,plot_data_fit_dfe,plot_data_fit_heatmap,data2mut_matrix,plot_data_comparison_bar,plot_data_fit_clustermap,plot_sub_matrix,plot_cov,plot_data_comparison_violin
from dms2dfe.lib.io_mut_files import getwildtypecov
from dms2dfe.lib.plot_pdb import vector2bfactor
from dms2dfe.lib.global_vars import mut_types_form

def main(prj_dh):
    """
    This module makes vizualization for output of `dms2dfe`.

    #. Scatter grid plots raw counts in replicates, if present.
    #. Mutation matrix. of frequencies of mutants (log scaled). 
    #. Scatter plots of raw counts among selected and unselected samples 
    #. Mutation matrix. of Fitness values. 
    #. DFE plot. ie. Distribution of Fitness values for samples.
    #. Projections on PDB. Average of fitness values per residue are projected onto PDB file.  

    :param prj_dh: path to project directory.
    """
    logging.info("start")
    if not exists(prj_dh) :
        logging.error("Could not find '%s'" % prj_dh)
        sys.exit()
    configure.main(prj_dh)
    from dms2dfe.tmp import info
    pdb_fh=info.pdb_fh
    norm_type=info.norm_type
    Ni_cutoff=int(info.Ni_cutoff)
    cctmr=info.cctmr
    
    if cctmr != 'nan':
        cctmr=[int("%s" % i) for i in cctmr.split(" ")]
        cctmr=[(cctmr[0],cctmr[1]),(cctmr[2],cctmr[3])]
    else:
        cctmr=None

    
    data_feats=pd.read_csv(prj_dh+"/data_feats/aas/feats_all")

    #tmp files
    plot_pdb_chimera_fhs_fh='%s/tmp/plot_pdb_chimera_fhs' % abspath(dirname(__file__))
    plot_pdb_chimera_fhs_f = open(plot_pdb_chimera_fhs_fh, 'w+')

    if not exists(prj_dh+"/plots"):
        makedirs(prj_dh+"/plots")

#plot for coverage
    plot_type="coverage"
    logging.info("processing: plot type: %s" % plot_type)
    lbls=pd.read_csv("%s/cfg/lbls" % prj_dh).set_index("varname",drop=True)
    type_form="aas"
    for lbl in lbls.index.values:
        plot_fh="%s/plots/%s/fig_%s_%s.pdf" % (prj_dh,type_form,plot_type,lbl)
        if not exists(plot_fh):
            try:
                data_cov=getwildtypecov(lbl,lbls,cctmr)
                data_lbl=pd.read_csv("%s/data_lbl/aas/%s" % (prj_dh,lbl))
                plot_cov(data_cov,data_lbl,plot_fh=plot_fh)
            except:
                logging.info("coverage not be plotted plot: sequencing data not provided.")

#plot data_lbl
    plot_type="repli"
    logging.info("processing: plot type: %s" % plot_type)
    plot_data_lbl_repli(prj_dh,lim_min=np.log2(Ni_cutoff))
    
    plot_type="heatmap"
    cbar_label="Frequency ($N_{i}$)"
    logging.info("processing: plot type: %s" % plot_type)    
    data_lbls=glob("%s/data_lbl/aas/*" % prj_dh)+glob("%s/data_lbl/cds/*" % prj_dh)
    data_lbls= [basename(fh) for fh in data_lbls]
    if len(data_lbls)!=0:
        for data_lbli in np.unique(data_lbls):
            for type_form in mut_types_form:
                if not exists("%s/plots/%s" % (prj_dh,type_form)):
                    makedirs("%s/plots/%s" % (prj_dh,type_form))
                data_lbl_fh = "%s/data_lbl/%s/%s" % (prj_dh,type_form,data_lbli)
                if exists(data_lbl_fh):
                    data_lbl=pd.read_csv(data_lbl_fh)
                    if "Unnamed: 0" in data_lbl.columns:
                        data_lbl=data_lbl.drop("Unnamed: 0", axis=1)
                    #heatmaps
                    plot_fh="%s/plots/%s/fig_%s_%s.pdf" % (prj_dh,type_form,plot_type,data_lbli)
                    if not exists(plot_fh):
                        plot_data_fit_heatmap(data_lbl,type_form,col='NiAcutlog',cmap="Blues",\
                                                        center=None,data_feats=data_feats,plot_fh=plot_fh,cbar_label=cbar_label)
                    else:
                        logging.info("already processed: %s" % basename(plot_fh))

#plot data_fit
    data_fits=glob("%s/data_fit/aas/*" % prj_dh)+glob("%s/data_fit/cds/*" % prj_dh)
    data_fits= [basename(fh) for fh in data_fits if "_WRT_" in basename(fh)]
    if len(data_fits)!=0:
        for data_fiti in np.unique(data_fits):
            for type_form in mut_types_form:
                if not exists("%s/plots/%s" % (prj_dh,type_form)):
                    makedirs("%s/plots/%s" % (prj_dh,type_form))

                data_fit_pair=  data_fiti.split('_WRT_')
                # print data_fit_pair

                data_fit_sel  =data_fit_pair[0]
                data_fit_unsel=data_fit_pair[1]
                data_fit_fh = "%s/data_fit/%s/%s" % (prj_dh,type_form,data_fiti)
                if exists(data_fit_fh):
                    data_fit=pd.read_csv(data_fit_fh)
                    if "Unnamed: 0" in data_fit.columns:
                        data_fit=data_fit.drop("Unnamed: 0", axis=1)

                    if len(data_fit.loc[:,"FiA"].unique())>10:
                        if  ("NiAunsel" in data_fit.columns.tolist()) and \
                            (len(data_fit.loc[:,"NiAunsel"].unique())>10) and \
                            (len(data_fit.loc[:,"NiAsel"].unique())>10):

                            logging.info("processing: %s/%s" % (type_form,data_fiti))

                            plot_type="scatter"
                            logging.info("processing: plot type: %s" % plot_type)    
                            plot_fh="%s/plots/%s/fig_%s_%s.pdf" % (prj_dh,type_form,plot_type,data_fiti) 
                            if not exists(plot_fh):
                                plot_data_fit_scatter(data_fit,norm_type,Ni_cutoff,plot_fh=plot_fh)
                            else:
                                logging.info("already processed: %s" % basename(plot_fh))

                        plot_type="dfe"
                        logging.info("processing: plot type: %s" % plot_type)    
                        plot_fh="%s/plots/%s/fig_%s_%s.pdf" % (prj_dh,type_form,plot_type,data_fiti) 
                        if not exists(plot_fh):
                            plot_data_fit_dfe(data_fit,norm_type,axvspan_min=0,axvspan_max=0,plot_fh=plot_fh)
                        else:
                            logging.info("already processed: %s" % basename(plot_fh))
                        
                        plot_type="heatmap"
                        logging.info("processing: plot type: %s" % plot_type)    
                        plot_fh="%s/plots/%s/fig_%s_%s.pdf" % (prj_dh,type_form,plot_type,data_fiti) 
                        if not exists(plot_fh):
                            if "_inferred" in basename(plot_fh):
                                cbar_label="Fitness scores ($F_{i}$)"
                            else:
                                cbar_label="Log fold change ($FC_{i}$)"
                            plot_data_fit_heatmap(data_fit,type_form,col='FiA',data_feats=data_feats,\
                                                  plot_fh=plot_fh,cbar_label=cbar_label)
                        else:
                            logging.info("already processed: %s" % basename(plot_fh))
                        
                        plot_type="clustermap_rows"
                        logging.info("processing: plot type: %s" % plot_type)    
                        plot_fh="%s/plots/%s/fig_%s_%s.pdf" % (prj_dh,type_form,plot_type,data_fiti) 
                        if not exists(plot_fh):
                            plot_data_fit_clustermap(data_fit,type_form,col='FiA',plot_fh=plot_fh)
                        else:
                            logging.info("already processed: %s" % basename(plot_fh))

                        # plot_type="clustermap_rows_cols"
                        # plot_fh="%s/plots/%s/fig_%s_%s.pdf" % (prj_dh,type_form,plot_type,data_fiti) 
                        # if not exists(plot_fh):
                        #     ax=plot_data_fit_clustermap(data_fit,type_form,col='FiA',col_cluster=True)
                        #     ax.savefig(plot_fh+".pdf",format='pdf')                        
                        #     ax.savefig(plot_fh);plt.clf();plt.close()
                        # else:
                        #     logging.info("already processed: %s" % basename(plot_fh))

                        plot_type="sub_matrix"
                        logging.info("processing: plot type: %s" % plot_type)    
                        plot_fh="%s/plots/%s/fig_%s_%s.pdf" % (prj_dh,type_form,plot_type,data_fiti) 
                        if not exists(plot_fh):
                            ax=plot_sub_matrix(data_fit,type_form,col='FiA',plot_fh=plot_fh)
                        else:
                            logging.info("already processed: %s" % basename(plot_fh))
                       
                        plot_type="pdb"
                        logging.info("processing: plot type: %s" % plot_type)    
                        if "aas" in type_form:
                            mut_matrix=data2mut_matrix(data_fit,'FiA','mut',type_form)
                            data_fit_avg=mut_matrix.mean()
                            pdb_clrd_fh="%s/plots/%s/fig_%s_%s.pdb" % (prj_dh,type_form,plot_type,data_fiti) 
                            if not exists(pdb_clrd_fh):
                                vector2bfactor(data_fit_avg,pdb_fh,pdb_clrd_fh)
                                plot_pdb_chimera_fhs_f.write(abspath(pdb_clrd_fh)+"\n")
                            else:
                                logging.info("already processed: %s" % basename(pdb_clrd_fh))
                    else:
                        logging.info("skipping data_fit : %s/%s : need more numbers" % (type_form,data_fiti))                   
                else:
                    logging.info("can not find : %s" % (data_fit_fh))
    else:
        logging.info("already processed")
    plot_pdb_chimera_fhs_f.close()

# data_comparison    
    data_comparisons=glob("%s/data_comparison/aas/*" % prj_dh)+glob("%s/data_comparison/cds/*" % prj_dh)
    data_comparisons= [basename(fh) for fh in data_comparisons]
    if len(data_comparisons)!=0:
        for data_comparisoni in np.unique(data_comparisons):
            for type_form in mut_types_form:
                if type_form=="aas":
                    if not exists("%s/plots/%s" % (prj_dh,type_form)):
                        makedirs("%s/plots/%s" % (prj_dh,type_form))
                    data_comparison_fh = "%s/data_comparison/%s/%s" % (prj_dh,type_form,data_comparisoni)
                    data_comparison=pd.read_csv(data_comparison_fh)
                    if "Unnamed: 0" in data_comparison.columns:
                        data_comparison=data_comparison.drop("Unnamed: 0", axis=1)
                    #bar   
                    plot_type="bar"
                    classes=["negative","robust","positive"]
                    plot_fh="%s/plots/%s/fig_%s_%s.pdf" % (prj_dh,type_form,plot_type,data_comparisoni) 
                    if not exists(plot_fh):
                        ax=plot_data_comparison_bar(data_comparison,plot_fh=plot_fh,index=classes)
                    else:
                        logging.info("already processed: %s" % basename(plot_fh))
                    plot_type="violin"
                    plot_fh="%s/plots/%s/fig_%s_%s.pdf" % (prj_dh,type_form,plot_type,data_comparisoni) 
                    if not exists(plot_fh):
                        plot_data_comparison_violin(data_comparison,plot_fh=plot_fh)
                    else:
                        logging.info("already processed: %s" % basename(plot_fh))

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

if __name__ == '__main__':
    main(sys.argv[1])