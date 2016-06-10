from os import makedirs
from os.path import exists
import pandas as pd
import subprocess
import logging
logging.basicConfig(format='[%(asctime)s] %(levelname)s\tfrom %(filename)s in %(funcName)s(..): %(message)s',level=logging.DEBUG) 
from dms2dfe.lib.plot_mut_data import data2mut_matrix

fn="Data_S1-S4.xlsx"
if not exists("TEM1_Firnberg_et_al_2014/data_input"):
    makedirs("TEM1_Firnberg_et_al_2014/data_input")
if not exists("TEM1_Firnberg_et_al_2014/data_input/Data_S1-S4.xlsx"):
    logging.info("Downloading supporting data")
    com="wget -q http://mbe.oxfordjournals.org/content/suppl/2014/02/22/msu081.DC1/"+fn+" --directory-prefix=TEM1_Firnberg_et_al_2014/data_input"
    subprocess.call(com,shell=True)
logging.info("Generating mutation-matrices (.mut_cds) in TEM1_Firnberg_et_al/data_input")
data=pd.read_excel("TEM1_Firnberg_et_al_2014/data_input/"+fn,"S1 Codon fitnesses")
cols_amp_concs=list(data.iloc[0,:])
data=data.iloc[1:,:]
cols=data.columns.tolist()
for i,col in enumerate(cols_amp_concs):
    if not pd.isnull(col):
        cols[i]=col
data.columns=cols
data.to_csv("TEM1_Firnberg_et_al_2014/data_input/"+fn+".csv",index=False)
cols_amp_concs=[col for col in cols_amp_concs if not pd.isnull(col)]
data.loc[:,"ref"]=data.loc[:,"WT codon"]
data.loc[:,"mut"]=data.loc[:,"Mutant codon"]
for i in data.index.values:
    data.loc[i,"mutids"]="%s%03d%s" % (data.loc[i,"ref"],data.loc[i,"Ambler Position"],data.loc[i,"mut"])

for col in [0.5,256,512]:#cols_amp_concs:
    if not exists("TEM1_Firnberg_et_al_2014/data_input/Amp%s.mat_mut_cds" % (col)):
        mut_matrix=data2mut_matrix(data, col, "mut", "cds")
        mut_matrix=mut_matrix.T
        mut_matrix=mut_matrix.reset_index()
        for i in range(len(mut_matrix)):
            mut_matrix.loc[i,"ref_cd"]= mut_matrix.loc[i,"refrefi"][:3]
        del mut_matrix["refrefi"]
        mut_matrix.to_csv("TEM1_Firnberg_et_al_2014/data_input/Amp%s.mat_mut_cds" % (col),index=False)
logging.info("Generated mut-matrix in TEM1_Firnberg_et_al/data_input/*.mut_cds")