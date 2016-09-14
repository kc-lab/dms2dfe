#!usr/bin/python

# Copyright 2016, Rohan Dandage <rraadd_8@hotmail.com,rohan@igib.in>
# This program is distributed under General Public License v. 3.  

import sys
import numpy as np
from os.path import exists,splitext
from os import makedirs,stat
from Bio import SeqIO
import pandas as pd
from multiprocessing import Pool
import logging
logging.basicConfig(format='[%(asctime)s] %(levelname)s\tfrom %(filename)s in %(funcName)s(..): %(message)s',level=logging.DEBUG) # filename=cfg_xls_fh+'.log'
from dms2dfe import configure
from dms2dfe.lib.io_mut_files import getusable_lbls_list,getusable_fits_list,mut_mat_cds2data_lbl,data_lbl2data_fit,repli2data_lbl_avg #,makemutids,mat_cds2mat_aas,getNS,collate_cctmr,repli2avg,class_fit,repli2data_lbl_avg
from dms2dfe.lib.global_vars import mut_types_form

def main(prj_dh):
    """
    This modules converts mutation matrices (.mat files produced in upstream ana1_sam2mutmat module) and calculates the fitness values for samples.
    The output data is saved in `data_fit` format as described in :ref:`io`.
    
    :param prj_dh: path to project directory.
    """
    logging.info("start")
    if not exists(prj_dh) :
        logging.error("Could not find '%s'" % prj_dh)
        sys.exit()
    configure.main(prj_dh)

    global prj_dh_global,host,norm_type,fsta_seqlen,cctmr_global,output_dh,prj_dh_global,lbls,Ni_cutoff

    from dms2dfe.tmp import info
    fsta_fh=info.fsta_fh
    cctmr=info.cctmr
    host=info.host
    cores=info.cores
    norm_type=info.norm_type
    Ni_cutoff=int(info.Ni_cutoff)
    # SET global variables

    lbls=pd.read_csv(prj_dh+'/cfg/lbls')
    lbls=lbls.set_index('varname')
    prj_dh_global=prj_dh
    
    if cctmr != 'nan':
        cctmr=[int("%s" % i) for i in cctmr.split(" ")]
        cctmr_global=[(cctmr[0],cctmr[1]),(cctmr[2],cctmr[3])]
    else:
        cctmr_global=None

    with open(fsta_fh,'r') as fsta_data:
        for fsta_record in SeqIO.parse(fsta_data, "fasta") :
            fsta_id=fsta_record.id
            fsta_seq=str(fsta_record.seq) 
            fsta_seqlen=len(fsta_seq)
            logging.info("ref name : '%s', length : '%d' " % (fsta_record.id, fsta_seqlen))

    lbls_list=getusable_lbls_list(prj_dh)
    if len(lbls_list)!=0:
        pool_mut_mat_cds2data_lbl=Pool(processes=int(cores)) 
        pool_mut_mat_cds2data_lbl.map(pooled_mut_mat_cds2data_lbl,lbls_list)
        pool_mut_mat_cds2data_lbl.close(); pool_mut_mat_cds2data_lbl.join()
        # pooled_mut_mat_cds2data_lbl(lbls_list[0])
    else:
        logging.info("already processed: mut_mat_cds2data_lbl")
        
    repli2data_lbl_avg(prj_dh)
    fits_pairs_list    =getusable_fits_list(prj_dh)    
    if len(fits_pairs_list)!=0:
        pool_data_lbl2data_fit=Pool(processes=int(cores)) 
        pool_data_lbl2data_fit.map(pooled_data_lbl2data_fit,fits_pairs_list)
        pool_data_lbl2data_fit.close(); pool_data_lbl2data_fit.join()
        # pooled_data_lbl2data_fit(fits_pairs_list[0])
    else:
        logging.info("already processed: data_lbl2data_fit")
    logging.shutdown()

def pooled_mut_mat_cds2data_lbl(lbls_list_tp):
    """
    This function converts mutation matrix (produced from .mat) file into data_lbl format.
    
    :param lbls_list_tp: tuple with name (lbl) of sample and path to file with codon level mutation matrix (.mat_mut_cds).
    """
    lbli=lbls_list_tp[0]
    lbl_mat_mut_cds_fh=lbls_list_tp[1]
    logging.info("processing : %s" % (lbli))
    mut_mat_cds2data_lbl(lbli,lbl_mat_mut_cds_fh,host,prj_dh_global,fsta_seqlen,cctmr_global,Ni_cutoff)


def pooled_data_lbl2data_fit(fits_list):
    """
    This estimates Fitness (data_fit) from mutation data (data_lbl) in selcted and unselected samples.  

    :param fits_list: tuple with name (lbl) of input and selected samples.
    """
    unsel_lbl=fits_list[0]
    sel_lbl=fits_list[1]
    logging.info("processing : %s and %s" % (unsel_lbl,sel_lbl))
    data_lbl2data_fit(unsel_lbl,sel_lbl,norm_type,prj_dh_global,cctmr_global,lbls)

if __name__ == '__main__':
    main(sys.argv[1])