#!usr/bin/python

# Copyright 2016, Rohan Dandage <rraadd_8@hotmail.com,rohan@igib.in>
# This program is distributed under General Public License v. 3.  

"""
================================
``get_protein_features``
================================
"""
import pandas as pd
from os.path import exists,abspath,dirname,basename
from Bio.PDB import DSSP,PDBParser
from Bio.PDB.Polypeptide import PPBuilder
import numpy as np
import subprocess
import re
import pysam

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from os.path import basename, splitext,exists
from dms2dfe.lib.io_seq_files import fasta_nts2prt

from Bio import SeqIO,Seq,SeqRecord
from Bio.Alphabet import IUPAC
import subprocess

import warnings
warnings.simplefilter(action = "ignore") # , category = PDBConstructionWarning
import logging
logging.basicConfig(format='[%(asctime)s] %(levelname)s\tfrom %(filename)s in %(funcName)s(..): %(message)s',level=logging.DEBUG) # filename=cfg_xls_fh+'.log'
from dms2dfe.lib.global_vars import aas_21_3letter,secstruc_lbls

def getdssp_data(pdb_fh,dssp_fh):
    """
    This uses DSSP to get structural information.

    :param pdb_fh: path to .pdb file
    :param dssp_fh: path to DSSP source.
    """
    dssp_data_fh="%s/../tmp/dssp"% (abspath(dirname(__file__)))
    dssp_com="%s -i %s -o %s" % (dssp_fh,pdb_fh,dssp_data_fh)
    subprocess.call(dssp_com,shell=True)

    start=False
    dssp_f=open(dssp_data_fh, 'r')
    dssp_fwf_f=open(dssp_data_fh+'.fwf', 'w')
    for line in dssp_f:
        if "#" in line:
            start=True
        if start is True:
            dssp_fwf_f.write(line)
    dssp_f.close()
    dssp_fwf_f.close()

    dssp_data=pd.read_fwf(dssp_data_fh+'.fwf',widths=[5,5,2,2,3,4,1,1,2,4,4,1,4,7,1,4,6,1,4,6,1,4,6,1,4,2,6,6,6,6,6,7,7,7])

    dssp_data.columns=["junk","aasi","chain","ref_dssp","Secondary structure","Helix formation in helix types 3 4 and 5","bend","Chirality",\
                       "beta bridge labels","First residue of beta bridge","Second residue of beta bridge","Sheet of beta bridge",\
                       "Solvent accessibility",\
                      "Offset from residue to the partner in N-H-->O H-bond (1)","junk","Energy (kcal/mol) for this N-H-->O H-bond (1)",\
                      "Offset from residue to the partner in O-->H-N H-bond (1)","junk","Energy (kcal/mol) for this O-->H-N H-bond (1)",\
                      "Offset from residue to the partner in N-H-->O H-bond (2)","junk","Energy (kcal/mol) for this N-H-->O H-bond (2)",\
                      "Offset from residue to the partner in O-->H-N H-bond (2)","junk","Energy (kcal/mol) for this O-->H-N H-bond (2)",\
                       "junk","cosine of the angle between C=O of residue and C=O of previous residue",\
                      "kappa bond/bend angle","alpha torsion/dihedral angle",\
                       "phi torsion angle","psi torsion angle",\
                      "X coordinates of CA atom","Y coordinates of CA atom","Z coordinates of CA atom"]

    dssp_data=dssp_data.loc[dssp_data.loc[:,"chain"]=='A',:] # only chain A
    del dssp_data["junk"]
    del dssp_data["chain"]
    del dssp_data["ref_dssp"]
    dssp_data=dssp_data.set_index("aasi")
    return dssp_data

def pdb2dfromactivesite(pdb_fh,active_sites=[]):
    """
    This calculates distances between each ligand atom or optionally provided amino acids (sources) and each residue in the protein.
    
    :param pdb_fh: path to .pdb file.
    :param active_sites: optional list of residue numbers as sources. 
    """
    junk_residues = ["HOH"," MG","CA"," NA","SO4","IOD","NA","CL","GOL","PO4"]
    pdb_parser=PDBParser()
    pdb_data=pdb_parser.get_structure("pdb_name",pdb_fh)
    model = pdb_data[0]
    chainA = model["A"] #only a chain
    residues   = list(chainA.get_residues())
    ligands_residue_objs=[]
    for residue in chainA:
        if not residue.get_resname() in junk_residues:
            if not residue.get_resname() in aas_21_3letter: #only aas 
                ligands_residue_objs.append(residue)
            elif residue.id[1] in active_sites:
                ligands_residue_objs.append(residue)
            
    dfromligands=pd.DataFrame()
    for ligandi in range(len(ligands_residue_objs)):
        ligand_residue_obj=ligands_residue_objs[ligandi]
        for ligand_atom_obj in ligand_residue_obj:
            for residue in chainA:
                if residue.get_resname() in aas_21_3letter: #only aas 
                    dfromligands.loc[residue.id[1],"ref_pdb"]=residue.get_resname()
                    if not ligand_residue_obj.get_resname() in aas_21_3letter:
                        dfromligands.loc[residue.id[1],"Distance from Ligand: %s (ATOM: %s)" % (ligand_residue_obj.get_resname(), \
                                            ligand_atom_obj.get_name())]=ligand_residue_obj[ligand_atom_obj.get_name()]-residue["CA"]
                    else:
                        dfromligands.loc[residue.id[1],"Distance from active site residue: %s (ATOM: %s)" % (ligand_residue_obj.get_resname(), \
                                            ligand_atom_obj.get_name())]=ligand_residue_obj[ligand_atom_obj.get_name()]-residue["CA"]

    dfromligands.index.name="aasi"
    del dfromligands["ref_pdb"]
    return dfromligands

def get_consrv_score(fsta_fh,host,clustalo_fh,rate4site_fh):
    from Bio.Blast import NCBIWWW
    from Bio.Blast import NCBIXML
    from skbio import TabularMSA, Protein
    from Bio.Alphabet import IUPAC

    blast_method='blastp'
    blast_db="swissprot"
    
    # make prt fasta
    fsta_prt_fh="%s_prt%s" % (splitext(fsta_fh)[0],splitext(fsta_fh)[1])
    if not exists(fsta_prt_fh):
        prt_seq=fasta_nts2prt(fsta_fh,host=host)

    for fsta_data in SeqIO.parse(fsta_prt_fh,'fasta'):
        ref_id=fsta_data.id
        prt_seq=fsta_data.seq
        break

    blast_fh="%s_blastp.xml" % (splitext(fsta_fh)[0])
    blast_fasta_fh="%s.fasta" % (splitext(blast_fh)[0])
    msa_fh=blast_fasta_fh.replace("blastp","blastp_msa")
    if not exists(msa_fh):
        if not exists(blast_fasta_fh):
            if not exists(blast_fh):
                # get blast
                blast_handle = NCBIWWW.qblast(blast_method, blast_db, sequence=prt_seq)
                blast_results = blast_handle.read()
                blast_f=open(blast_fh, 'w')
                blast_f.write(blast_results)
                blast_f.close()
            blast_f = open(blast_fh, 'r')
            blast_records=NCBIXML.parse(blast_f)

            # blast to fasta
            blast_fasta_f=open(blast_fasta_fh, 'w')
            fsta_data=[]
            fsta_data.append(SeqRecord.SeqRecord(Seq.Seq(prt_seq,IUPAC.protein),id = ref_id,description=''))
            for rec in blast_records:
                for aln in rec.alignments:
                    for hsp in aln.hsps:
                        fsta_data.append(SeqRecord.SeqRecord(Seq.Seq(hsp.sbjct,IUPAC.protein),id = aln.hit_id,description=''))
            SeqIO.write(fsta_data, blast_fasta_f, "fasta")
            blast_fasta_f.close()
    # blast fasta to msa : clustaw
        clustalo_com="./%s -i %s -o %s --outfmt=fasta --log=%s.log" % (clustalo_fh,blast_fasta_fh,msa_fh,msa_fh)
        subprocess.call(clustalo_com,shell=True)

    # msa to consrv skbio==0.4.1
    msa = TabularMSA.read(msa_fh, constructor=Protein)
    msa.reassign_index(minter='id')
    metric='inverse_shannon_uncertainty'
    gap_modes=['ignore','include']
    gap_modes_labels=['gaps ignored','gaps included'] #'no gaps',
    data_feats_conserv=pd.DataFrame(columns=["aasi"])

    for fsta_data in SeqIO.parse(msa_fh,'fasta'):
        if ref_id in fsta_data.id:
            ref_seq=fsta_data.seq
            break

    for gap_mode in gap_modes:
        positional_conservation = msa.conservation(metric=metric, degenerate_mode='nan', gap_mode=gap_mode)
        aai=1
        for i in range(len(ref_seq)):
            if ref_seq[i]!="-":
                data_feats_conserv.loc[i,"aasi"]=aai
                data_feats_conserv.loc[i,"Conservation score (inverse shannon uncertainty): %s" % (gap_modes_labels[gap_modes.index(gap_mode)])]=positional_conservation[i]
                aai+=1
    data_feats_conserv=data_feats_conserv.set_index("aasi")
    
    rate4site_out_fh="%s.rate4site" % (splitext(rate4site_fh)[0])
    rate4site_com="./%s -s %s -o %s -a %s" % (rate4site_fh,msa_fh,rate4site_out_fh,ref_id)
    subprocess.call(rate4site_com,shell=True)
    rate4site_out_csv_fh="%s.rate4site.csv" % (splitext(rate4site_out_fh)[0])
    with open(rate4site_out_fh,"r") as rate4site_out_f:
        lines = rate4site_out_f.readlines()
    with open(rate4site_out_csv_fh,'w') as rate4site_out_csv_f:
        for line in lines:
            if (12<lines.index(line)<len(lines)-2):
                rate4site_out_csv_f.write(line[:19]+"\n")
    data_feats_conserv_rate4site=pd.read_csv(rate4site_out_csv_fh,delim_whitespace=True, header=None,names=["aasi","ref","Conservation score (rate4site)"])
    data_feats_conserv_rate4site=data_feats_conserv_rate4site.drop("ref",axis=1).set_index("aasi")
    
    data_feats_conserv=pd.concat([data_feats_conserv,data_feats_conserv_rate4site],axis=1)
    return data_feats_conserv


def get_residue_depth(pdb_fh,msms_fh):
    from Bio.PDB import Selection,PDBParser
    from Bio.PDB.Polypeptide import is_aa
    from Bio.PDB.ResidueDepth import get_surface,_read_vertex_array,residue_depth,ca_depth,min_dist
    surface_fh="%s/%s.msms.vert" % (dirname(msms_fh),basename(pdb_fh))
    if not exists(surface_fh):
        pdb_to_xyzr_fh="%s/pdb_to_xyzr" % dirname(msms_fh)
        xyzr_fh="%s/%s.xyzr" % (dirname(msms_fh),basename(pdb_fh))
        pdb_to_xyzr_com="./%s %s > %s" % (pdb_to_xyzr_fh,pdb_fh,xyzr_fh)
        msms_com="./%s -probe_radius 1.5 -if %s -of %s > %s.log" % (msms_fh,xyzr_fh,splitext(surface_fh)[0],splitext(surface_fh)[0])
        subprocess.call("%s;%s" % (pdb_to_xyzr_com,msms_com) , shell=True)
    surface =_read_vertex_array(surface_fh)
    
    pdb_parser=PDBParser()
    pdb_data=pdb_parser.get_structure("pdb_name",pdb_fh)
    model = pdb_data[0]
    residue_list = Selection.unfold_entities(model, 'R') 
    
    depth_dict = {}
    depth_list = []
    depth_keys = []
    for residue in residue_list:
        if not is_aa(residue):
            continue
        rd = residue_depth(residue, surface)
        ca_rd = ca_depth(residue, surface)
        # Get the key
        res_id = residue.get_id()
        chain_id = residue.get_parent().get_id()
        depth_dict[(chain_id, res_id)] = (rd, ca_rd)
        depth_list.append((residue, (rd, ca_rd)))
        depth_keys.append((chain_id, res_id))
        # Update xtra information
        residue.xtra['EXP_RD'] = rd
        residue.xtra['EXP_RD_CA'] = ca_rd

    depth_df=pd.DataFrame(depth_dict).T.reset_index()
    depth_df=depth_df.drop("level_0",axis=1)
    aasi_prev=0
    for i in range(len(depth_df)):
        if depth_df.loc[i,"level_1"][1]!=aasi_prev:
            depth_df.loc[i,"aasi"]=depth_df.loc[i,"level_1"][1]
            aasi_prev=depth_df.loc[i,"level_1"][1]

    depth_df=depth_df.drop("level_1",axis=1)
    depth_df=depth_df.loc[~pd.isnull(depth_df.loc[:,"aasi"]),:]
    depth_df=depth_df.set_index("aasi",drop=True)
    depth_df.columns=["Residue depth","Residue (C-alpha) depth"]
    return depth_df