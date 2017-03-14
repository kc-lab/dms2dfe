#!usr/bin/python

# Copyright 2016, Rohan Dandage <rraadd_8@hotmail.com,rohan@igib.in>
# This program is distributed under General Public License v. 3.  

import sys
from os import makedirs,stat
from os.path import splitext, join, exists, isdir,basename,abspath,dirname
import pandas as pd 

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
    if not exists(prj_dh+"/plots"):
        makedirs(prj_dh+"/plots")
        
    plot_coverage(info)
    plot_heatmap(info)
    plot_clustermap(info)
    plot_multisca(info)
    plot_pdb(info)
    plot_violin(info)
    plot_pies(info)
    
if __name__ == '__main__':
    main(sys.argv[1])