#!usr/bin/python

# Copyright 2016, Rohan Dandage <rraadd_8@hotmail.com,rohan@igib.in>
# This program is distributed under General Public License v. 3.  

import sys
from os import makedirs,stat
from os.path import splitext, join, exists, isdir,basename,abspath,dirname
import pandas as pd 
from dms2dfe.lib.io_plot_files import plot_coverage,plot_heatmap,plot_clustermap,plot_multisca,plot_violin,plot_pies,plot_pdb

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
    # fsta_fh=info.fsta_fh
    # pdb_fh=info.pdb_fh
    # norm_type=info.norm_type
    # Ni_cutoff=int(info.Ni_cutoff)
    # cctmr=info.cctmr            
    # if cctmr != 'nan':
    #     cctmr=[int("%s" % i) for i in cctmr.split(" ")]
    #     cctmr=[(cctmr[0],cctmr[1]),(cctmr[2],cctmr[3])]
    # else:
    #     cctmr=None
    # data_feats=pd.read_csv(prj_dh+"/data_feats/aas/feats_all")
    # #tmp files
    # plot_pdb_chimera_fhs_fh='%s/tmp/plot_pdb_chimera_fhs' % abspath(dirname(__file__))
    # plot_pdb_chimera_fhs_f = open(plot_pdb_chimera_fhs_fh, 'w+')
    if not exists(prj_dh+"/plots"):
        makedirs(prj_dh+"/plots")
        
    plot_coverage(info,plot_type="coverage")
    plot_heatmap(info,data_fit_fhs,plot_type="heatmap")
    plot_clustermap
    plot_multisca
    plot_violin
    plot_pies
    plot_pdb

if __name__ == '__main__':
    main(sys.argv[1])