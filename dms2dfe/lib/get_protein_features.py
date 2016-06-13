#!usr/bin/python

# Copyright 2016, Rohan Dandage <rraadd_8@hotmail.com,rohan@igib.in>
# This program is distributed under General Public License v. 3.  

"""
================================
``get_protein_features``
================================
"""
import pandas as pd
from os.path import exists,abspath,dirname
from Bio.PDB import DSSP,PDBParser
from Bio.PDB.Polypeptide import PPBuilder
import numpy as np
import subprocess
import re
import pysam
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