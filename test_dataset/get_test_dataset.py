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

import pandas as pd
import numpy as np
from scipy import stats
from os.path import exists,dirname,basename,isfile,splitext
from os import makedirs
import matplotlib.pyplot as plt
from glob import glob
import seaborn as sns
from dms2dfe.lib.io_nums import str2num
import logging
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)
import subprocess

from dms2dfe.lib.io_nums import is_numeric
from dms2dfe import configure,ana2_mutmat2fit,ana4_modeller
from dms2dfe.lib.io_mut_files import rescale_fitnessbysynonymous
from dms2dfe.lib.io_ml import denanrows
from dms2dfe.lib.plot_mut_data import data2mut_matrix

def get_frequencies(prj_dh):
    if not exists(prj_dh):
        makedirs(prj_dh)
    if not exists("%s/data_lbl/aas" % prj_dh):
        makedirs("%s/data_lbl/aas" % prj_dh)

    configure.main(prj_dh)

    #get frequencies of mutants in data_lbl format

    if "TEM1" in prj_dh:
        lbls=["Amp0.5","Amp256"]
        ref="TEM1"
    elif "APH2" in prj_dh:
        lbls=["KKA2_Bkg1","KKA2_S1_Kan14_L1"]
        ref="APH2"
    unselected=lbls[0]    
    subprocess.call("cp %s/cfg/feats %s/cfg/feats" % (prj_dh.replace("_random_sampling",""),prj_dh),shell=True)
    subprocess.call("cp %s/cfg/info %s/cfg/info" % (prj_dh.replace("_random_sampling",""),prj_dh),shell=True)
    subprocess.call("cp %s/cfg/fit %s/cfg/fit" % (prj_dh.replace("_random_sampling",""),prj_dh),shell=True)
    subprocess.call("cp -r %s/data_feats %s/data_feats" % (prj_dh.replace("_random_sampling",""),prj_dh),shell=True)
    fit=pd.read_csv("%s/cfg/fit" % (prj_dh))
    for col in fit.columns:
        for i in fit.index.values:
            if not is_numeric(fit.loc[i,col]):
                if "_from_seq_data" in fit.loc[i,col]:
                    fit.loc[i,col]=np.nan
    fit.to_csv("%s/cfg/fit" % (prj_dh),index=False)    
    info=pd.read_csv("%s/cfg/info" % (prj_dh)).set_index("varname",drop=True)
    info.loc["cores","input"]=6
    info.to_csv("%s/cfg/info" % (prj_dh))    
    for lbl in lbls:
        subprocess.call("cp %s/data_lbl/aas/%s %s/data_lbl/aas/%s" % (prj_dh.replace("_random_sampling",""),lbl,prj_dh,lbl),shell=True)
    return ref,unselected

def make_random_splits(prj_dh):
    fit=pd.read_csv("%s/cfg/fit" % (prj_dh))
    for i in range(len(fit)):
        data_unseln=fit.loc[i,"unsel"]
        data_seln=fit.loc[i,"sel_1"]
        rand_states=[1,2,3,4,5]
        for rand_state in rand_states:
            data_unsel_fh="%s/data_lbl/aas/%s" % (prj_dh,data_unseln)
            data_sel_fh="%s/data_lbl/aas/%s" % (prj_dh,data_seln)
            if exists(data_unsel_fh) and exists(data_sel_fh):
                if not "_random_sampling" in basename(data_unsel_fh):
                    data_unsel=pd.read_csv(data_unsel_fh).set_index("mutids",drop=True)
                    data_sel=pd.read_csv(data_sel_fh).set_index("mutids",drop=True)
                    if 'Unnamed: 0' in data_unsel.columns.tolist():
                        data_unsel=data_unsel.drop("Unnamed: 0",axis=1)
                    if 'Unnamed: 0' in data_sel.columns.tolist():
                        data_sel=data_sel.drop("Unnamed: 0",axis=1)

                    np.random.seed(rand_state)
                    random_rowis=np.random.randint(1,len(data_unsel),1000)
                    random_mutids=[]
                    for random_rowi in random_rowis:
                        random_mutid=data_sel.index.values[random_rowi]
                        if (not len(random_mutids)>=100) \
                        and (random_mutid not in random_mutids):
                            if not pd.isnull(data_sel.loc[random_mutid,"NiA"])\
                            and not pd.isnull(data_unsel.loc[random_mutid,"NiA"]):
                                random_mutids.append(random_mutid)
                    data_unsel_test=data_unsel.loc[random_mutids,:]
                    data_sel_test=data_sel.loc[random_mutids,:]

                    data_unsel_train=data_unsel.loc[:,:]
                    data_sel_train=data_sel.loc[:,:]
                    data_unsel_train.loc[random_mutids,["NiA","NiAcut","NiAcutlog","NiScutlog"]]=np.nan
                    data_sel_train.loc[random_mutids,["NiA","NiAcut","NiAcutlog","NiScutlog"]]=np.nan

                    if not exists("%s/random_sampling_tests" % prj_dh):
                        makedirs("%s/random_sampling_tests" % prj_dh)
                    data_unsel_test.reset_index().to_csv("%s/random_sampling_tests/%s_random_sampling_%d_test" % (prj_dh,data_unseln,rand_state),index=False)
                    data_sel_test.reset_index().to_csv("%s/random_sampling_tests/%s_random_sampling_%d_test" % (prj_dh,data_seln,rand_state),index=False)
                    data_unsel_train.reset_index().to_csv("%s/data_lbl/aas/%s_random_sampling_%d_train" % (prj_dh,data_unseln,rand_state),index=False)
                    data_sel_train.reset_index().to_csv("%s/data_lbl/aas/%s_random_sampling_%d_train" % (prj_dh,data_seln,rand_state),index=False)    
                    fit.loc[len(fit)+i,"unsel"]="%s_random_sampling_%d_train" % (data_unseln,rand_state)
                    fit.loc[len(fit)+i-1,"sel_1"]="%s_random_sampling_%d_train" % (data_seln,rand_state)
    fit.to_csv("%s/cfg/fit" % (prj_dh),index=False)

def make_dmstools_inputs(prj_dh):
    if not exists("%s/data_dms_tools" % prj_dh):
        makedirs("%s/data_dms_tools" % prj_dh)

    for fh in glob("%s/data_lbl/aas/*" % prj_dh):
        if isfile(fh):
            data=pd.read_csv(fh)
            col="NiA"
            mut_matrix=data2mut_matrix(data, col, "mut", "aas")
            mut_matrix=mut_matrix.T
            mut_matrix=mut_matrix.reset_index()
            for i in range(len(mut_matrix)):
                mut_matrix.loc[i,"WT"]= mut_matrix.loc[i,"refrefi"].replace("X","*")[:1]
                mut_matrix.loc[i,"# POSITION"]= str2num(mut_matrix.loc[i,"refrefi"])
            del mut_matrix["refrefi"]
            mut_matrix=mut_matrix.loc[:,["# POSITION","WT","A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","*"]]
            mut_matrix_fh="%s/data_dms_tools/%s.dms_tools.input" % (prj_dh,basename(fh))
            mut_matrix.fillna(0).to_csv(mut_matrix_fh,sep=" ",index=False)

            with open(mut_matrix_fh, 'r') as f:
                line = f.read()
                line=line.replace("\"","")
            with open(mut_matrix_fh, 'w') as f:
                f.write(line)

def run_dmstools(prj_dh,dms_tools_fh,unselected):
    input_fhs=glob("%s/data_dms_tools/*.dms_tools.input" % (prj_dh))
    selected_fhs=[fh for fh in input_fhs if not unselected in fh]
    unselected_fhs=[fh for fh in input_fhs if unselected in fh]

    unsel_sel_tuples=zip(np.sort(unselected_fhs),np.sort(selected_fhs))
    for i in range(len(unsel_sel_tuples)):
        unseli=unsel_sel_tuples[i][0]
        seli  =unsel_sel_tuples[i][1]
    #     com="python %s %s %s %s_WRT_%s.dms_tools.prefs --ncpus 8 --chartype aa &\n" % (dms_tools_fh,unseli,seli,seli,splitext(basename(unseli))[0])
        com="python %s %s %s %s_WRT_%s.dms_tools.prefs --chartype aa &\n" % (dms_tools_fh,unseli,seli,seli,splitext(basename(unseli))[0])
        print com

def get_dmstools_output(prj_dh,ref):
    inferprefs_fhs=glob("%s/data_dms_tools/*.dms_tools.prefs" % (prj_dh))

    for inferprefs_fh in inferprefs_fhs:
        if not exists(inferprefs_fh+".dms2dfe.data_fit"):
            with open(inferprefs_fh, 'r') as f:
                line = f.read()
                line=line.replace("# POSI","POSI")
                line=line.replace("*","X")
            with open(inferprefs_fh+"_mod", 'w') as f:
                f.write(line)

            inferprefs=pd.read_csv(inferprefs_fh+"_mod",sep=" ")
            mut_matrix=inferprefs.loc[:,['POSITION',
             'WT',
             'PI_A',
             'PI_C',
             'PI_D',
             'PI_E',
             'PI_F',
             'PI_G',
             'PI_H',
             'PI_I',
             'PI_K',
             'PI_L',
             'PI_M',
             'PI_N',
             'PI_P',
             'PI_Q',
             'PI_R',
             'PI_S',
             'PI_T',
             'PI_V',
             'PI_W',
             'PI_Y',
             'PI_X' ]]

            mut_matrix.columns=[col.replace("PI_","") for col in mut_matrix.columns.tolist()]
            mut_matrix=mut_matrix.sort_values(by="POSITION",axis=0)
            mut_matrix["ref"]=mut_matrix["WT"]
            mut_matrix=mut_matrix.set_index("ref")
            mut_matrix.columns.name='mut'
            mut_matrix.index.name='ref'

            data_fit=pd.DataFrame(mut_matrix.loc[:,['A',
             'C',
             'D',
             'E',
             'F',
             'G',
             'H',
             'I',
             'K',
             'L',
             'M',
             'N',
             'P',
             'Q',
             'R',
             'S',
             'T',
             'V',
             'W',
             'Y',
             'X']].unstack())

            data_fit.columns=["pref"]
            data_fit=data_fit.reset_index()
            data_fit_dms_tools=data_fit

            # data_fit_dms_tools.to_csv("/home/kclabws1/Documents/propro/writ/prjs/1_dms_software/data/dms_tools/BTEM1_ana/TEM1_Firnberg_et_al_2014/data_fit/Amp256_WRT_Amp0.5",index=False)

            # add ref i
            ref_prev="B"

            for i in range(len(data_fit_dms_tools)):
                if data_fit_dms_tools.loc[i,"mut"] != ref_prev:
                    ref_prev=data_fit_dms_tools.loc[i,"mut"]
                    if ref=="TEM1":
                        refi=3
                    else:
                        refi=1
                data_fit_dms_tools.loc[i,"refi"] = refi    
                data_fit_dms_tools.loc[i,"refrefi"] = "%s%03d" % (data_fit_dms_tools.loc[i,"ref"],refi) 
                data_fit_dms_tools.loc[i,"mutids"] = "%s%03d%s" % (data_fit_dms_tools.loc[i,"ref"],refi,data_fit_dms_tools.loc[i,"mut"])

                if ref=="TEM1":
                    if refi==238 or refi==252:
                        refi+=2
                    else:
                        refi+=1
                else:
                    refi+=1

            data_fit_dms_tools.loc[:,"FCA"]=data_fit_dms_tools.loc[:,"pref"]
            data_fit_dms_tools=rescale_fitnessbysynonymous(data_fit_dms_tools)

            data_fit_dms_tools.loc[data_fit_dms_tools.loc[:,"FiA"]>1,"class_fit"]="beneficial"
            data_fit_dms_tools.loc[data_fit_dms_tools.loc[:,"FiA"]<1,"class_fit"]="deleterious"
            data_fit_dms_tools.loc[data_fit_dms_tools.loc[:,"FiA"]==1,"class_fit"]="neutral"

            data_fit_dms_tools.to_csv(inferprefs_fh+".dms2dfe.data_fit",index=False)

def get_comparison_data(prj_dh):
    fit=pd.read_csv("%s/cfg/fit" % (prj_dh))
    for i in fit.index.values:
        if not pd.isnull(fit.loc[i,"unsel"]):
            if "train" in fit.loc[i,"unsel"]:
                data_unsel_trainn=fit.loc[i,"unsel"]
                data_sel_trainn  =fit.loc[i,"sel_1"]
                data_unseln=data_unsel_trainn[:data_unsel_trainn.index("_rando")]
                data_seln  =data_sel_trainn[:data_sel_trainn.index("_rando")]
                data_random_sampling_compare_fh="%s/random_sampling_tests/%s_WRT_%s.compare" % (prj_dh,data_sel_trainn,data_unsel_trainn)
                if not exists(data_random_sampling_compare_fh):
                    data_unsel_test =pd.read_csv("%s/random_sampling_tests/%s" % (prj_dh,data_unsel_trainn.replace("train","test"))).set_index("mutids",drop=True)
                    data_fit=pd.read_csv("%s/data_fit/aas/%s_WRT_%s" % (prj_dh,data_seln,data_unseln)).set_index("mutids",drop=True)

                    data_random_sampling_compare=pd.DataFrame()

                    data_fit_train_infered_d2=pd.read_csv("%s/data_fit/aas/%s_WRT_%s_inferred" % (prj_dh,data_sel_trainn,data_unsel_trainn)).set_index("mutids",drop=True)
                    data_random_sampling_compare.loc[:,"Fold changes (true)"]=data_fit.loc[data_unsel_test.index.values,"FCA"]
                    data_random_sampling_compare.loc[:,"$FC_{i}$ (dms2dfe)"]=data_fit_train_infered_d2.loc[data_unsel_test.index.values,"FCA"]            
                    data_random_sampling_compare.loc[:,"$F_{i}$ (dms2dfe)"]=data_fit_train_infered_d2.loc[data_unsel_test.index.values,"FiA"]            

                    data_fit_train_infered_dt=pd.read_csv("%s/data_dms_tools/%s.dms_tools.input_WRT_%s.dms_tools.dms_tools.prefs.dms2dfe.data_fit" % (prj_dh,data_sel_trainn,data_unsel_trainn)).set_index("mutids",drop=True)
                    data_random_sampling_compare.loc[:,"Fold changes (true)2"]=data_fit.loc[data_unsel_test.index.values,"FCA"]
                    data_random_sampling_compare.loc[:,"$\pi_{r,x}$ (dms_tools)"]=data_fit_train_infered_dt.loc[data_unsel_test.index.values,"FCA"]
                    data_random_sampling_compare.loc[:,"$\Phi_{r,x}$ (dms_tools)"]=data_fit_train_infered_dt.loc[data_unsel_test.index.values,"FiA"]

                    data_random_sampling_compare.to_csv(data_random_sampling_compare_fh)

def make_comparison_data(prj_dh):
    fit=pd.read_csv("%s/cfg/fit" % (prj_dh))
    fit=denanrows(fit,condi="all").reset_index()
    plt.figure(figsize=[12,5],dpi=300)
    ax1=plt.subplot(121)
    ax2=plt.subplot(122)
    colors=['gray','r','g','b','m','c','k']
    legends1=[]
    legends2=[]
    for i in fit.index.values:
        if "train" in fit.loc[i,"unsel"]:
            data_unsel_trainn=fit.loc[i,"unsel"]
            data_sel_trainn  =fit.loc[i,"sel_1"]
            data_unseln=data_unsel_trainn[:data_unsel_trainn.index("_rando")]
            data_seln  =data_sel_trainn[:data_sel_trainn.index("_rando")]
            data_random_sampling_compare_fh="%s/random_sampling_tests/%s_WRT_%s.compare" % (prj_dh,data_sel_trainn,data_unsel_trainn)

            data_random_sampling_compare=pd.read_csv(data_random_sampling_compare_fh)
            sns.regplot(data=data_random_sampling_compare,x="Fold changes (true)",y="$F_{i}$ (dms2dfe)",color=colors[i],ax=ax1)
            data_random_sampling_compare=denanrows(data_random_sampling_compare)
            r, _ = stats.pearsonr(data_random_sampling_compare.loc[:,"Fold changes (true)"],\
                                 data_random_sampling_compare.loc[:,"$F_{i}$ (dms2dfe)"])
            legends1.append("random sample %d: r=%.2f" % (len(legends1)+1,r))

            data_random_sampling_compare=pd.read_csv(data_random_sampling_compare_fh)
            sns.regplot(data=data_random_sampling_compare,x="Fold changes (true)",y="$\Phi_{r,x}$ (dms_tools)",color=colors[i],ax=ax2)
            data_random_sampling_compare=denanrows(data_random_sampling_compare)
            r, _ = stats.pearsonr(data_random_sampling_compare.loc[:,"Fold changes (true)"],\
                                 data_random_sampling_compare.loc[:,"$\Phi_{r,x}$ (dms_tools)"])
            legends2.append("random sample %d: r=%.2f" % (len(legends2)+1,r))

    ax1.set_xlabel("$FC_{i}$ (empirical)")
    ax2.set_xlabel("$FC_{i}$ (empirical)")
    ax1.legend(legends1,loc="lower right")
    # ax2.legend(legends2,loc="upper right")
    ax2.legend(legends2,loc="lower right")
    plot_fh="%s/fig_corr_combo.pdf" % dirname(data_random_sampling_compare_fh)
    plt.savefig(plot_fh,format='pdf')
    print "Output plot: %" % plot_fh
    
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
    data_input_fns=["SRX547509/SRR1292901_1.fastq","SRX547509/SRR1292901_2.fastq",
                    "SRX547496/SRR1292881_1.fastq","SRX547496/SRR1292881_2.fastq",
                    "SRX547413/SRR1292709_1.fastq","SRX547413/SRR1292709_2.fastq"]   
    if not exists(data_input_dh):
        makedirs(data_input_dh)
    log_fh="%s.get_test_dataset.log" % data_input_dh
    log_f = open(log_fh,'a')
    for data_input_fn in data_input_fns:
        data_input_fh="%s/%s" % (data_input_dh,data_input_fn)
        if not exists(data_input_fh):
            logging.info("downloading: %d/6 files" % (data_input_fns.index(data_input_fn)+1))
            com="wget -q ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/fastq/SRA082/SRA082074/%s.bz2 --directory-prefix=%s" % (data_input_fn,data_input_dh)
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
        data_type=raw_input("Get 'publication data' or 'sequencing data (961.3 MB)'? [valid inputs: 1 | 2]:")
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