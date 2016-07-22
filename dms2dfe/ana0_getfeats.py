#!usr/bin/python

# Copyright 2016, Rohan Dandage <rraadd_8@hotmail.com,rohan@igib.in>
# This program is distributed under General Public License v. 3.  

import sys
from os.path import exists,splitext
from os import makedirs,stat
import logging
import numpy as np
import pandas as pd
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
from Bio import SeqIO
import logging
logging.basicConfig(format='[%(asctime)s] %(levelname)s\tfrom %(filename)s in %(funcName)s(..): %(message)s',level=logging.DEBUG) # filename=cfg_xls_fh+'.log'
from dms2dfe import configure
from dms2dfe.lib.io_seq_files import cctmr_fasta2ref_fasta
from dms2dfe.lib.io_data_files import convert2h5form
from dms2dfe.lib.get_protein_features import getdssp_data,pdb2dfromactivesite,get_consrv_score,get_residue_depth

def main(prj_dh):
    """
    This optional module extracts following structural features of reference protein.
    
    The out files are created in `prj_dh/data_feats`

    The steps and required dependendencies are following. 
    
    .. code-block:: text
    
        Secondary structure                      : using DSSP.
        Solvent Accessible Surface Area          : using DSSP.  
        Distance of a residue from reference atom: using Bio.PDB

    :param prj_dh: path to project directory.
    """
    logging.basicConfig(format='[%(asctime)s] %(levelname)s from %(funcName)s:\t%(message)s',level=logging.DEBUG) # filename=cfg_xls_fh+'.log'
    logging.info("start")

    if not exists(prj_dh) :
        logging.error("Could not find '%s'\n" % prj_dh)
        sys.exit()
    configure.main(prj_dh)
    from dms2dfe.tmp import info
    cctmr=info.cctmr
    fsta_fh=info.fsta_fh
    pdb_fh=info.pdb_fh
    dssp_fh=info.dssp_fh
    active_sites=info.active_sites
    host=info.host
    clustalo_fh=info.clustalo_fh
    msms_fh=info.msms_fh
    rate4site_fh=info.rate4site_fh
    
    if active_sites!='nan':
        active_sites=[int(i) for i in active_sites.split(" ")]
    else:
        active_sites=[]

    if cctmr != 'nan':
        cctmr=[int("%s" % i) for i in cctmr.split(" ")]
        aas_len=cctmr[1]-1
        fsta_fh=cctmr_fasta2ref_fasta(fsta_fh,cctmr)
    else :
        fsta_data = SeqIO.read(open(fsta_fh), "fasta")
        aas_len=len(fsta_data)/3

    if not exists(prj_dh+"/data_feats/aas"):
        makedirs(prj_dh+"/data_feats/aas")

    if not exists("%s/data_feats/feats" % prj_dh):
        if exists("%s/cfg/feats" % prj_dh) and stat("%s/cfg/feats" % prj_dh).st_size !=0:
            data_feats=pd.read_csv('%s/cfg/feats' % prj_dh)
            if len(data_feats)!=0:
                if "Unnamed: 0" in data_feats.columns:
                    data_feats=data_feats.drop("Unnamed: 0", axis=1)
                data_feats=data_feats.set_index("aasi",drop=True)
                data_feats.index = [int(i) for i in data_feats.index]
                tmp=pd.DataFrame(index=range(1,aas_len+1,1))
                tmp.index.name="aasi"
                data_feats=pd.concat([tmp,data_feats],axis=1,join_axes=[tmp.index])
                logging.info("custom features taken from: cfg/feats")
        else:
            data_feats=pd.DataFrame(index=range(1,aas_len+1,1))
            data_feats.index.name="aasi"
        if not pd.isnull(pdb_fh):
            feats_types=["dssp","dfromact","depth","consrv_score"]
            for feats_type in feats_types:
                data_feats_fh="%s/data_feats/aas/feats_%s" % (prj_dh,feats_type)
                if feats_type=="dssp":
                    if not exists(data_feats_fh):
                        dssp_df=getdssp_data(pdb_fh,dssp_fh)
                        dssp_df.reset_index().to_csv(data_feats_fh,index=False)
                    else:
                        dssp_df=pd.read_csv(data_feats_fh)
                        dssp_df=dssp_df.set_index("aasi",drop=True)
                elif feats_type=="dfromact":
                    if not exists(data_feats_fh):
                        dfromact_df=pdb2dfromactivesite(pdb_fh,active_sites)
                        dfromact_df.reset_index().to_csv(data_feats_fh,index=False)
                    else:
                        dfromact_df=pd.read_csv(data_feats_fh)
                        dfromact_df=dfromact_df.set_index("aasi",drop=True)                        
                elif feats_type=="depth":
                    logging.info("getting structural features")
                    if not exists(data_feats_fh):
                        depth_df=get_residue_depth(pdb_fh,msms_fh)
                        depth_df.reset_index().to_csv(data_feats_fh,index=False)
                    else:
                        depth_df=pd.read_csv(data_feats_fh)
                        depth_df=depth_df.set_index("aasi",drop=True)                        
                elif feats_type=="consrv_score":
                    logging.info("getting conservation scores")
                    if not exists(data_feats_fh):
                        consrv_score_df=get_consrv_score(fsta_fh,host,clustalo_fh,rate4site_fh)
                        consrv_score_df.reset_index().to_csv(data_feats_fh,index=False)            
                    else:
                        consrv_score_df=pd.read_csv(data_feats_fh)
                        consrv_score_df=consrv_score_df.set_index("aasi",drop=True)                        
            if len(data_feats)!=0:
                data_feats=pd.concat([data_feats,dssp_df,dfromact_df,consrv_score_df,depth_df], axis=1) #combine dssp_df and dfromact
            else:
                data_feats=pd.concat([dssp_df,dfromact_df,consrv_score_df,depth_df], axis=1) #combine dssp_df and dfromact
            data_feats.reset_index().to_csv("%s/data_feats/aas/feats_all" % prj_dh,index=False)
            logging.info("output: data_feats/aas/feats_all")
        else:
            logging.warning("pdb_fh not given in cfg")
    else:
        logging.info("already processed")
    logging.shutdown()

if __name__ == '__main__':
    main(sys.argv[1])