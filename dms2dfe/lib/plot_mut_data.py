#!usr/bin/python

# Copyright 2016, Rohan Dandage <rraadd_8@hotmail.com,rohan@igib.in>
# This program is distributed under General Public License v. 3.  

"""
================================
``plot``
================================
"""
import sys
from os.path import splitext,exists,basename
from os import makedirs,stat
import pandas as pd
import numpy as np
import matplotlib 
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
# matplotlib.style.use('ggplot')
import logging
logging.basicConfig(format='[%(asctime)s] %(levelname)s\tfrom %(filename)s in %(funcName)s(..): %(message)s',level=logging.DEBUG) # 
from dms2dfe.lib.global_vars import mut_types_form


def plot_data_lbl_repli(prj_dh):
    import seaborn as sns
    """
    This plots scatters of frequecies of replicates.
    
    :param prj_dh: path to project directory.   
    """
    repli=pd.read_csv('%s/cfg/repli' % prj_dh)
    if "Unnamed: 0" in repli.columns:
        repli=repli.drop("Unnamed: 0", axis=1)

    for i,replii in repli.iterrows():
        replii=[str(lbl) for lbl in replii if not 'nan' in str(lbl)] # if not 'nan' in lbl
        if len(replii)!=0:
            for type_form in mut_types_form:
                if not exists("%s/plots/%s" % (prj_dh,type_form)):
                    makedirs("%s/plots/%s" % (prj_dh,type_form))
                data_lbl_key="data_lbl/%s/%s" % (type_form,replii[0])
                plot_fh="%s/plots/%s/fig_repli_%s.pdf" % (prj_dh,type_form,data_lbl_key.replace('/','_'))
                if not exists(plot_fh): 
                    data_repli=pd.DataFrame()
                    data_repli_usable=[]
                    for lbl in replii:
                        data_lbl_fh= "%s/data_lbl/%s/%s" % (prj_dh,type_form,lbl)
                        try:
                            data_lbl=pd.read_csv(data_lbl_fh)
                            data_repli[lbl]=data_lbl['NiAcutlog'].fillna(0)
                            if len(data_repli.loc[:,lbl].unique())>10:
                                data_repli_usable.append(True)
                            else:
                                data_repli_usable.append(False)
                        except:
                            logging.warning("do not exist: %s" % lbl)
                    data_repli=data_repli.loc[:,data_repli_usable]
                    if len(data_repli.columns)!=0:
                        ax = sns.pairplot(data_repli, kind="reg")
                        plt.tight_layout()
                        ax.savefig(plot_fh,format="pdf")
                        plt.clf();plt.close()
                        logging.info("output: %s" % basename(plot_fh))
                    else:
                        logging.warning("skipping data_replii : %s" % (data_lbl_key))
                else:
                    logging.info("already processed: %s" % basename(plot_fh))

def plot_data_fit_scatter(data_fit,norm_type,Ni_cutoff,plot_fh=None):
    """
    This plots scatter plot of frequecies of selected and unselected samples.
    
    :param data_fit: fitness data (`data_fit`). 
    :param norm_type: Type of normalization across samples [wild: wrt wild type | syn : wrt synonymous mutations | none : fold change serves as fitness]
    :param Ni_cutoff: threshold of depth per codon to be processed further. 
    """
    if np.isinf(np.log2(Ni_cutoff)):
        ax_lim_min=0 
    else:
        ax_lim_min=np.log2(Ni_cutoff)
    ax_lim_max=np.ceil(np.max([data_fit.loc[:,'NiAunsel'].max(),data_fit.loc[:,'NiAsel'].min()*-1]))+3
    
    fig = plt.figure(figsize=(3, 3),dpi=500)
    ax=plt.subplot(111)
    # plt.gcf().subplots_adjust(bottom=0.15)
    if norm_type=="syn":
        data_fit.plot(kind='scatter', x='NiAunsel', y='NiAsel', \
                           color='Black', label='Non-synonymous', ax=ax)
        data_fit.plot(kind='scatter', x='NiSunsel', y='NiSsel', \
                      color='darkorange', label='Synonymous', ax=ax)
    else:
        data_fit.plot(kind='scatter', x='NiAunsel', y='NiAsel', \
                           color='Black', label=None, ax=ax)    
    ax.set_xlabel(r'$N_{i,unsel}$')
    ax.set_ylabel(r'$N_{i,sel}$')
    ax.set_xlim([ax_lim_min,ax_lim_max])
    ax.set_ylim([ax_lim_min,ax_lim_max])
    ax.legend(loc= "upper right",frameon=True)
    plt.tight_layout()
    if plot_fh!=None:
        plt.savefig(plot_fh+".pdf",format='pdf')                        
    return ax

def plot_data_fit_dfe(data_fit,norm_type,col_fit="FiA",xlabel=r'$F_{i}$',plot_fh=None):
    """
    This plots histogram of Distribution of Fitness Effects (DFE).
    
    :param data_fit: fitness data (`data_fit`).
    :param norm_type: Type of normalization across samples [wild: wrt wild type | syn : wrt synonymous mutations | none : fold change serves as fitness]
    """
    fig = plt.figure(figsize=(8, 3),dpi=500)
    ax=plt.subplot(111)
    if norm_type=="syn":
        if not 'FiS' in data_fit.columns.tolist():
                data_fit.loc[:,'FiS']=data_fit.loc[(data_fit_infered.loc[:,'ref']==data_fit.loc[:,'mut']),col_fit]
        
        data_fit[[col_fit,'FiS']].plot(kind='hist',bins=40,color=['seagreen','darkorange'],\
                                     ax=ax)
        l2=ax.legend(['Non-synonymous','Synonymous'],loc="upper right")
    else:
        data_fit[col_fit].plot(kind='hist',bins=40,color='seagreen',ax=ax)

    ax.axvspan(-20, -2, color='blue', alpha=0.10)
    ax.axvspan(-2, 2, color='grey', alpha=0.10)
    ax.axvspan(2, 20, color='red', alpha=0.10)
    # l1=ax.legend(['Deleterious','Neutral','Beneficial'], loc='upper left')

    ax.set_xlabel(xlabel)
    ax.set_ylabel('Count')
    xlim=np.ceil(np.max([data_fit.loc[:,col_fit].max(),data_fit.loc[:,col_fit].min()*-1]))
    ax.set_xlim(-xlim,xlim)
    plt.tight_layout()
    if plot_fh!=None:
        plt.savefig(plot_fh+".pdf",format='pdf')                        
    return ax

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
    data_fit_heatmap=pd.pivot_table(data_fit,values=values_col,index=index_col,columns='refrefi')
    data_fit_heatmap=data_fit_heatmap.reindex_axis(data_fit.loc[:,'refrefi'].unique(),axis='columns')
    return data_fit_heatmap

def plotsspatches(ssi,ssi_ini,aai,aai_ini,patches,patches_colors,ax):
# H Alpha helix
# B Beta bridge
# E Strand
# G Helix-3
# I Helix-5
# T Turn
# S Bend
    sss=["H","B","E","G","I","T","S"]
    for ss in sss:
        if ssi==ss:
            ssi_ini=ssi
            aai_ini=aai+1
        elif ssi_ini==ss:
            width=aai+1-aai_ini
            if ss == "H":
                patch=mpatches.Arrow(aai_ini, 0, width,0, 4, edgecolor="none")
            else:
                patch=mpatches.Rectangle([aai_ini, -0.5], width, 1,edgecolor="none")                
            patches.append(patch)
            patches_colors.append(sss.index(ss))
            ssi_ini=ssi
            aai_ini=aai
            ax.text(aai_ini-width/2,0.5,ss,fontdict={'size': 16})
    return patches,patches_colors,aai_ini,ssi_ini,ax

def plotss(data_feats,ax):
#     plt.figure(figsize=(40, 1),dpi=500)
#     ax=plt.subplot(111)
    ax.set_axis_off()
    xlim=len(data_feats)

    #no sec struct
    for aai,rowi in data_feats.iterrows():
        if not pd.isnull(rowi.loc["Secondary structure"]):
            ini=aai+1
            break
    for aai,rowi in data_feats.iloc[::-1].iterrows():
        if not pd.isnull(rowi.loc["Secondary structure"]):
            end=aai+1
            break
    x,y = np.array([[ini, end], [0, 0]])
    line = mlines.Line2D(x, y, lw=5,color="black")
    line.set_zorder(0)
    ax.add_line(line)

    patches = []
    patches_colors=[]
    ssi_ini=""
    aai_ini=0
    for aai,rowi in data_feats.fillna("").iterrows():
        aai=aai-1
        ssi=rowi.loc["Secondary structure"]
        if ssi!=ssi_ini:
            patches,patches_colors,aai_ini,ssi_ini,ax=plotsspatches(ssi,ssi_ini,aai,aai_ini,patches,patches_colors,ax)
    collection = PatchCollection(patches, cmap=plt.cm.Accent, alpha=1)
    colors = np.linspace(0, 1, len(patches))
    collection.set_array(np.array(patches_colors))
    collection.set_zorder(20)
    ax.add_collection(collection)
    # ax.set_xticklabels(range(xlim))
    ax.set_xlim((0-0.5,xlim+0.5))
    ax.set_ylim((-1.1,1.1))
    ax.text(xlim+1,-0.5,"Secondary structure",fontdict={'size': 20})
#     plt.show()
    return ax

def plotacc(data_feats,ax):
    ax.set_axis_off()
    xlim=len(data_feats)
    ylim=data_feats.loc[:,"Solvent accessibility"].max()
    ax.stackplot(data_feats.loc[:,"aasi"]-1,data_feats.loc[:,"Solvent accessibility"])
    ax.set_xlim((0-0.5,xlim+0.5))
    ax.set_ylim((0,ylim))
    ax.text(xlim+1,0,"Solvent accessibility",fontdict={'size': 20})
    return ax    

def plot_data_fit_heatmap(data_fit,type_form,col,cmap="coolwarm",center=0,data_feats=None,xticklabels=None,plot_fh=None):
    """
    This plots heatmap of fitness values.
    
    :param data_fit: input data (`data_lbl` or `data_fit`) (dataframe).
    :param type_form: type of values ["aas" : amino acid | "cds" : codons].
    :param col: eg. columns with values. col for data_fit
    :param cmap: name of colormap (str).
    :param center: center colormap to this value (int).
    :param data_feats: input features `data_feats`.
    :param xticklabels: xticklabels of heatmap ["None" : reference index | "seq" : reference sequence].
    """
    from dms2dfe.lib.io_nums import str2num
    from dms2dfe.lib.global_vars import aas_21,cds_64
    import seaborn as sns
    
    data_fit_heatmap  =data2mut_matrix(data_fit,col,'mut',type_form)

    refis=[str2num(i) for i in data_fit_heatmap.columns.tolist()]
    refrefis=pd.DataFrame(data_fit_heatmap.columns.tolist(),index=refis,columns=["refrefi"])

    data_fit_heatmap2=pd.DataFrame(index=data_fit_heatmap.index)
    for i in range(1,refis[-1]+1):
        if i not in refis:
            data_fit_heatmap2.loc[:,i]=np.nan
        else :
            data_fit_heatmap2.loc[:,refrefis.loc[i,"refrefi"]]=data_fit_heatmap.loc[:,refrefis.loc[i,"refrefi"]]

    data_syn_locs=data_fit.loc[0:len(data_fit)/21-1,["mutids","ref"]]
    data_syn_locs["refi"]=[str2num(i)-1+0.15 for i in data_syn_locs["mutids"]]
    if "aas" in type_form:
        data_syn_locs["muti"]=[20-aas_21.index(i)+0.15 for i in data_syn_locs["ref"]]
    if "cds" in type_form:
        cds_64.sort()
        data_syn_locs["muti"]=[63-cds_64.index(i)+0.15 for i in data_syn_locs["ref"]]


    fig=plt.figure(figsize=(80, 12),dpi=500)      
    gs = gridspec.GridSpec(3, 1,height_ratios=[1,1,32])


    ax_all=plt.subplot(gs[:])
    ax_all.set_axis_off()

    ax = plt.subplot(gs[2])
    result=sns.heatmap(data_fit_heatmap2,cmap=cmap,ax=ax)
    ax.set_xlabel('Wild type',fontdict={'size': 20})
    ax.set_ylabel('Mutation to',fontdict={'size': 20})
    cbar=ax.figure.colorbar(ax.collections[0])
    cbar.set_label(("$%s$" % col),fontdict={'size': 20})
    if xticklabels=="seq":
        ax.set_xticklabels(data_fit_heatmap2.columns.tolist(),rotation=90)
    else:
        ax.set_xticks(range(1,len(data_fit_heatmap2.columns),1))
        ax.set_xticklabels(range(1,len(data_fit_heatmap2.columns),1),rotation=90)
    yticklabels=data_fit_heatmap2.index.values.tolist()
    ax.set_yticklabels(yticklabels[::-1],rotation=0)

    for i in data_syn_locs.index.values:
        ax.text(data_syn_locs.loc[i,"refi"],data_syn_locs.loc[i,"muti"],"s")

    if not data_feats is None: 
        ax_ss = plt.subplot(gs[0])
        ax_acc = plt.subplot(gs[1])

        ax_pos=ax.get_position()
        ax_ss_pos=ax_ss.get_position()
        ax_ss.set_position([ax_ss_pos.x0,ax_ss_pos.y0-0.05,ax_pos.width,ax_ss_pos.height*2])
        # ax_ss=plt.axes([ax_ss_pos.x0,ax_ss_pos.y0,ax_pos.width,ax_ss_pos.height])
        ax_ss.set_axis_off()
        ax_acc_pos=ax_acc.get_position()
        ax_acc.set_position([ax_acc_pos.x0,ax_acc_pos.y0-0.03,ax_pos.width,ax_acc_pos.height*2])
        ax_acc.set_axis_off()

        ax_ss=plotss(data_feats,ax_ss)
        ax_acc=plotacc(data_feats,ax_acc)
    extent = ax_all.get_window_extent().transformed(ax_all.figure.dpi_scale_trans.inverted())
    extent.set_points(np.array([[5,0],[36,11]]))    
    
    if plot_fh!=None:
        ax_all.figure.savefig(plot_fh+".pdf",format='pdf', bbox_inches=extent)                        
        # ax_all.figure.savefig(plot_fh, bbox_inches=extent);plt.clf();plt.close()
    return ax_all,extent


def plot_data_comparison_bar(data_comparison):
    """
    This plots the proportion of mutants according to type of selection in effect.
    
    :param data_comparison: input `data_comparison` (dataframe).
    """
    fig = plt.figure(figsize=(3, 3))
    ax=plt.subplot(111)

    index=["positive","negative","robust"]

    plot_data=pd.DataFrame(index=index,columns=["count"])
    group_data=pd.DataFrame(data_comparison.groupby("class_comparison").count().loc[:,"mutids"])
    for i in index:
        if i in group_data.index.values:
            plot_data.loc[i,"count"]=group_data.loc[i,"mutids"]

    plot_data.plot(kind="bar",ax=ax)
    ax.set_xlabel('Classes of comparison')
    ax.set_ylabel('Count')
    ax.legend().set_visible(False)
    plt.tight_layout()
    return ax

def plot_data_fit_clustermap(data_fit,type_form,col,cmap="coolwarm",center=0,col_cluster=False,row_cluster=True,plot_fh=None):
    """
    This clusters heatmaps.
    
    :param data_fit: fitness data
    :param type_form: type of mutants ["aas" : amino acid | "cds" : codon level]
    :param col: eg. columns with values. col for data_fit
    """
    import seaborn as sns
    data_fit_heatmap  =data2mut_matrix(data_fit,col,'mut',type_form)
    plt.figure()
    ax=sns.clustermap(data_fit_heatmap.fillna(0),method='average', metric='euclidean',\
                      col_cluster=col_cluster,row_cluster=row_cluster)
    ax.ax_heatmap.set_xlabel('Wild type')
    ax.ax_heatmap.set_ylabel('Mutation to')
    if not col_cluster:
        ax.ax_heatmap.set_xticks(range(1,len(data_fit_heatmap.columns),20))
        ax.ax_heatmap.set_xticklabels(range(1,len(data_fit_heatmap.columns),20),rotation=90)
    if plot_fh!=None:
        plt.savefig(plot_fh+".pdf",format='pdf')                        
    return ax

def data2sub_matrix(data_fit,values_col,index_col,type_form): 
    """
    This creates substitition matrix from input data (frequncies `data_lbl` or `data_fit`).
    
    :param data_fit: fitness data (`data_fit`).
    :param values_col: column name with values (str). 
    :param index_col: column name with index (str).
    :param type_form: type of values ["aas" : amino acid | "cds" : codons].
    """                  
    sub_matrix=pd.pivot_table(data_fit,values=values_col,index=index_col,columns='ref')
    return sub_matrix

def plot_sub_matrix(data_fit,type_form,col,cmap="coolwarm",center=0,plot_fh=None):
    """
    This plots heatmap of fitness values.
    
    :param data_fit: input data (`data_lbl` or `data_fit`) (dataframe).
    :param type_form: type of values ["aas" : amino acid | "cds" : codons].
    :param col: eg. columns with values. col for data_fit
    :param cmap: name of colormap (str).
    :param center: center colormap to this value (int).
    """
    from dms2dfe.lib.io_nums import str2num
    from dms2dfe.lib.global_vars import aas_21,cds_64
    import seaborn as sns
    
    sub_matrix  =data2sub_matrix(data_fit,col,'mut',type_form)
    plt.figure(figsize=(5,4),dpi=500)
    ax=plt.subplot(111)
    result=sns.heatmap(sub_matrix.T,cmap=cmap,ax=ax)
    ax.set_ylabel('Wild type')
    ax.set_xlabel('Mutation to')
    yticklabels=sub_matrix.index.values.tolist()
    ax.set_yticklabels(yticklabels[::-1],rotation=0)
    plt.tight_layout()
    if plot_fh!=None:
        plt.savefig(plot_fh+".pdf",format='pdf')                        
    return ax


def plot_cov(data_cov,data_lbl,plot_fh=None):
    plt.figure(figsize=[6,3],dpi=300)
    ax1=plt.subplot(111)
    ax1.plot(data_cov,lw=2,color='b')
    ax2 = ax1.twinx()
    ax2.plot(data2mut_matrix(data_lbl,"NiA","mut","aas").sum().tolist(),lw=2,color='r')
    ax1.set_xlabel("Codon position")
    ax1.set_ylabel("Coverage (reads per codons)",color='b')
    ax2.set_ylabel("Mutants per codons",color='r')
    ax1.set_xlim([data_cov.index.values[0],data_cov.index.values[-1]])
    ax1.set_yscale('log')
    ax2.set_yscale('log')
    for tl in ax1.get_yticklabels():
        tl.set_color('b')
    for tl in ax2.get_yticklabels():
        tl.set_color('r')
    plt.tight_layout()
    if plot_fh!=None:
        plt.savefig(plot_fh+".pdf",format='pdf')                        
    return ax1,ax2