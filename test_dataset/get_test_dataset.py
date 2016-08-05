import numpy as np
from os.path import basename,dirname,splitext,exists
from os import stat
from os import makedirs
import pandas as pd
import subprocess
import logging
logging.basicConfig(format='[%(asctime)s] %(levelname)s\tfrom %(filename)s in %(funcName)s(..): %(message)s',level=logging.DEBUG) 
from glob import glob
from dms2dfe.lib.io_nums import str2num 

def get_TEM1_dataset(prj_dh):
    def data2mut_matrix(data_fit,values_col,index_col,type_form): 
        """
        This creates mutation matrix from input data (frequncies `data_lbl` or `data_fit`).

        :param data_fit: fitness data (`data_fit`).
        :param values_col: column name with values (str). 
        :param index_col: column name with index (str).
        :param type_form: type of values ["aas" : amino acid | "cds" : codons].
        """
        if 'aas' in type_form:   
            data_fit.loc[:,'refrefi']=(data_fit.mutids.astype(str).str[0:4])                                    
        elif 'cds' in type_form:                
            data_fit.loc[:,'refrefi']=(data_fit.mutids.astype(str).str[0:6])                    
        data_fit_heatmap=pd.pivot_table(data_fit,values=values_col,index=index_col,columns='refrefi',aggfunc='sum')
        data_fit_heatmap=data_fit_heatmap.reindex_axis(data_fit.loc[:,'refrefi'].unique(),axis='columns')
        return data_fit_heatmap
    fn="Data_S1-S4.xlsx"
    data_input_dh="%s/data_input" % prj_dh
    data_input_fh="%s/%s" % (data_input_dh,fn)      
    if not exists(data_input_dh):
        makedirs(data_input_dh)
    if not exists(data_input_fh):
        log_fh="%s.log" % data_input_fh
        log_f = open(log_fh,'a')
        com="wget -q http://mbe.oxfordjournals.org/content/suppl/2014/02/22/msu081.DC1/%s --directory-prefix=%s" % (fn,data_input_dh)
        subprocess.call(com,shell=True,stdout=log_f, stderr=subprocess.STDOUT)
        log_f.close()

        if not exists("%s/data_input/%s.csv" % (prj_dh,fn)):
            data=pd.read_excel(data_input_fh,"S1 Codon fitnesses")
            cols_amp_concs=list(data.iloc[0,:])
            data=data.iloc[1:,:]
            cols=data.columns.tolist()
            for i,col in enumerate(cols_amp_concs):
                if not pd.isnull(col):
                    cols[i]=col
            data.columns=cols
            data.to_csv("%s/data_input/%s.csv" % (prj_dh,fn),index=False)
            cols_amp_concs=[col for col in cols_amp_concs if not pd.isnull(col)]
            data.loc[:,"ref"]=data.loc[:,"WT AA"]
            data.loc[:,"mut"]=data.loc[:,"Mutant AA"]

            for i in data.index.values:
                data.loc[i,"mutids"]="%s%03d%s" % (data.loc[i,"ref"],data.loc[i,"Ambler Position"],data.loc[i,"mut"])
            for col in [0.5,256]:#cols_amp_concs:
                if not exists("%s/data_lbl/aas/Amp%s" % (prj_dh,col)):
                    data_aacounts=data.loc[:,["mutids",col]].reset_index(drop=True)
                    data_aacounts.columns=["mutids","NiA"]
                    for i in range(len(data_aacounts)):
                        if "*" in data_aacounts.loc[i,"mutids"]:
                            data_aacounts.loc[i,"mutids"]=data_aacounts.loc[i,"mutids"].replace("*","X")
                        data_aacounts.loc[i,"mut"]="%s" % (data_aacounts.loc[i,"mutids"][-1])
                    data_aacounts=data_aacounts.fillna(0)
                    mut_matrix=data2mut_matrix(data_aacounts, "NiA", "mut", "aas")
                    data_aacounts=mut_matrix.unstack().reset_index()
                    data_aacounts.columns=["refrefi","mut","NiA"]
                    data_aacounts.loc[data_aacounts.loc[:,"NiA"]==0,"NiA"]=np.nan
                    for i in range(len(data_aacounts)):
                        data_aacounts.loc[i,"ref"]="%s" % (data_aacounts.loc[i,"refrefi"][0])
                        data_aacounts.loc[i,"mutids"]="%s%s" % (data_aacounts.loc[i,"refrefi"],data_aacounts.loc[i,"mut"])
                        data_aacounts.loc[i,"NiAcut"]=data_aacounts.loc[i,"NiA"]
                        data_aacounts.loc[i,"NiAcutlog"]=np.log2(data_aacounts.loc[i,"NiA"])
                        if data_aacounts.loc[i,"mut"]==data_aacounts.loc[i,"ref"]:
                            data_aacounts.loc[i,"NiScutlog"]=data_aacounts.loc[i,"NiAcutlog"] 
                    if not exists("%s/data_lbl/aas" % (prj_dh)):
                        makedirs("%s/data_lbl/aas" % (prj_dh))
                    data_aacounts.to_csv("%s/data_lbl/aas/Amp%s" % (prj_dh,col))

def get_APH2_dataset(prj_dh):
    data_input_dh="%s/data_input" % prj_dh
    data_input_fh="%s/nar-01049-met-h-2014-File007.zip" % data_input_dh
    if not exists(data_input_dh):
        makedirs(data_input_dh)
    if not exists(data_input_fh):
        log_fh="%s.log" % data_input_fh
        log_f = open(log_fh,'a')
        com="wget -q http://nar.oxfordjournals.org/content/suppl/2014/06/09/gku511.DC1/nar-01049-met-h-2014-File007.zip --directory-prefix=%s" % (data_input_dh)
        subprocess.call(com,shell=True,stdout=log_f, stderr=subprocess.STDOUT)
        subprocess.call("unzip %s -d %s" % (data_input_fh,data_input_dh),shell=True,stdout=log_f, stderr=subprocess.STDOUT)
        log_f.close()

    if not exists("%s/data_lbl/aas" % prj_dh):
        makedirs("%s/data_lbl/aas" % prj_dh)
    aacounts_fns=["KKA2_Bkg1.aacounts.txt",
                    "KKA2_S1_Kan14_L1.aacounts.txt"]
    for fn in aacounts_fns:
        fh="%s/data_input/Supplementary Data 2 - Unpacked/%s" % (prj_dh,fn)
        data_lbl_fh="%s/data_lbl/aas/%s" % (prj_dh,splitext(splitext(basename(fh))[0])[0])
        if not exists(data_lbl_fh):
            data_aacounts=pd.read_csv(fh,header=1,sep="\t").set_index("Position",drop=True).T.reset_index()
            for i in range(len(data_aacounts)):
                data_aacounts.loc[i,"refrefi"]="%s%03d" % (data_aacounts.loc[i,"Wild-type"],str2num(data_aacounts.loc[i,"index"]))
            data_aacounts=data_aacounts.loc[:,[ 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', 'X', 'refrefi']]
            data_aacounts=data_aacounts.set_index('refrefi',drop=True)
            data_aacounts=data_aacounts.unstack().reset_index()
            data_aacounts.columns=["mut","refrefi","NiA"]
            for i in range(len(data_aacounts)):
                data_aacounts.loc[i,"ref"]="%s" % (data_aacounts.loc[i,"refrefi"][0])
                data_aacounts.loc[i,"mutids"]="%s%s" % (data_aacounts.loc[i,"refrefi"],data_aacounts.loc[i,"mut"])
                if not pd.isnull(data_aacounts.loc[i,"NiA"]):
                    data_aacounts.loc[i,"NiA"]=str2num(data_aacounts.loc[i,"NiA"])
            data_aacounts.loc[data_aacounts.loc[:,"NiA"]==0,"NiA"]=np.nan
            for i in range(len(data_aacounts)):
                data_aacounts.loc[i,"NiAcut"]=data_aacounts.loc[i,"NiA"]
                data_aacounts.loc[i,"NiAcutlog"]=np.log2(data_aacounts.loc[i,"NiA"])
                if data_aacounts.loc[i,"mut"]==data_aacounts.loc[i,"ref"]:
                    data_aacounts.loc[i,"NiScutlog"]=data_aacounts.loc[i,"NiAcutlog"]    
            data_aacounts.to_csv(data_lbl_fh)

def get_APH2_seq_data(prj_dh):
    data_input_dh="%s/data_input" % prj_dh
    data_input_fns=["SRR1292901_1.fastq","SRR1292901_2.fastq",
                    "SRR1292881_1.fastq","SRR1292881_2.fastq",
                    "SRR1292709_1.fastq","SRR1292709_2.fastq"]   
    if not exists(data_input_dh):
        makedirs(data_input_dh)
    log_fh="%s.get_test_dataset.log" % data_input_dh
    log_f = open(log_fh,'a')
    for data_input_fn in data_input_fns:
        data_input_fh="%s/%s" % (data_input_dh,data_input_fn)
        if not exists(data_input_fh):
            logging.info("downloading: %d/6 files" % (data_input_fns.index(data_input_fn)+1))
            com="wget -q ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/fastq/SRA082/SRA082074/SRX547509/%s.bz2 --directory-prefix=%s" % (data_input_fn,data_input_dh)
            # print com
            subprocess.call(com,shell=True,stdout=log_f, stderr=subprocess.STDOUT)
    subprocess.call("bzip2 -d %s/*.bz2" % data_input_dh ,shell=True,stdout=log_f, stderr=subprocess.STDOUT)
    log_f.close()
    
def main():
    test_dataset=raw_input("Get input data for TEM1 or APH2 dataset? [valid inputs: 1 | 2]:")
    if (test_dataset=='1') or ("TEM1" in test_dataset):
        logging.info("Getting 'TEM1' dataset..")
        prj_dh="TEM1_Firnberg_et_al_2014"
        get_TEM1_dataset(prj_dh)
    elif (test_dataset=='2') or ("APH2" in test_dataset):
        logging.info("Getting 'APH2' dataset..")
        data_type=raw_input("Get 'publication data' or 'sequencing data'? [valid inputs: 1 | 2]:")
        prj_dh="APH2_Melnikov_et_al_2014"    
        if data_type=='1':
            get_APH2_dataset(prj_dh)
        elif data_type=='2':    
            get_APH2_seq_data(prj_dh)
        else:
            logging.warning("[valid inputs: 1 | 2]")
    else:
        logging.warning("[valid inputs: 1 | 2]")

if __name__ == '__main__':
    main()