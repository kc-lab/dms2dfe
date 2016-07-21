#!usr/bin/python

# Copyright 2016, Rohan Dandage <rraadd_8@hotmail.com,rohan@igib.in>
# This program is distributed under General Public License v. 3.  

"""
================================
``io_ml``
================================
"""
from os.path import abspath,dirname,exists

from sklearn.ensemble import RandomForestClassifier
from sklearn.cross_validation import train_test_split
from sklearn.preprocessing import LabelEncoder,label_binarize
from sklearn.metrics import roc_curve, auc
from sklearn.grid_search import GridSearchCV

from dms2dfe.lib.io_nums import is_numeric

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import warnings
warnings.simplefilter(action = "ignore", category = FutureWarning)
import logging
logging.basicConfig(format='[%(asctime)s] %(levelname)s\tfrom %(filename)s in %(funcName)s(..): %(message)s',level=logging.DEBUG) # filename=cfg_xls_fh+'.log'

def data_fit_feats2combo(data_fit,data_feats,y_coln,keep_mutids=False):
    """
    This combines data_fit and data_feats to make data_all.
    
    :param data_fit: fitness data (pandas dataframe).
    :param data_feats: features wrt length of protein 
    :param y_coln: name of column with classes in `data_fit`. 
    :returns data_all: `data_fit` and `data_feats` concatenated.
    :returns Xcols: list of features.
    """

    if "Unnamed: 0" in data_feats.columns:
        data_feats=data_feats.drop("Unnamed: 0", axis=1)
    if "Unnamed: 0" in data_fit.columns:
        data_fit=data_fit.drop("Unnamed: 0", axis=1)

    for row in data_fit.iterrows():
        mutid=row[1]["mutids"]
        data_fit.loc[row[0],"aasi"]=int(''.join(ele for ele in mutid if ele.isdigit()))
    data_fit.loc[:,"aasi"]=data_fit.loc[:,"aasi"].astype(int)
    data_fit=data_fit.set_index("aasi",drop=True)
    
    if not data_feats.index.name=="aasi":
        data_feats=data_feats.set_index("aasi",drop=True)
    if np.nan in data_feats.index.values:
        data_feats=data_feats.drop(np.nan,axis=0)
    # data_feats.to_csv("test_data_feats")
    data_feats_prt=pd.DataFrame(index=data_fit.index,columns=data_feats.columns)
    rowi=0
    for row in data_feats_prt.iterrows():
        if (row[0] in data_feats_prt.index.values) and (row[0] in data_feats.index.values):
            data_feats_prt.iloc[rowi,:]=data_feats.loc[row[0],:]
        rowi+=1
    # data_feats.to_csv("test_data_feats2")
    data_feats_prt.columns=["%s" % col for col in data_feats_prt.columns.tolist()]
    # data_feats_prt.to_csv("test_data_feats_prt")
    
    data_feats_aas_fh='%s/data_feats_aas' % abspath(dirname(__file__))
    data_feats_aas=pd.read_csv(data_feats_aas_fh)
    data_feats_aas=data_feats_aas.set_index("aas",drop=True)

    data_feats_aas_mut=pd.DataFrame(index=data_fit.loc[:,"mut"],columns=data_feats_aas.columns)
    rowi=0
    for row in data_feats_aas_mut.iterrows():
        if row[0] in data_feats_aas.index.values:
            data_feats_aas_mut.iloc[rowi,:]=data_feats_aas.loc[row[0],:]
        rowi+=1
    data_feats_aas_mut.columns=["Mutant amino acid's %s" % col for col in data_feats_aas_mut.columns.tolist()]

    data_feats_aas_ref=pd.DataFrame(index=data_fit.loc[:,"ref"],columns=data_feats_aas.columns)
    rowi=0
    for row in data_feats_aas_ref.iterrows():
        if row[0] in data_feats_aas.index.values:
            data_feats_aas_ref.iloc[rowi,:]=data_feats_aas.loc[row[0],:]
        rowi+=1
    data_feats_aas_ref.columns=["Reference amino acid's %s" % col for col in data_feats_aas_ref.columns.tolist()]

    if len(data_fit) == len(data_feats_prt):
        if len(data_fit) == len(data_feats_aas_mut):
            if len(data_fit) == len(data_feats_aas_ref):
                data_feats_prt.index=data_fit.index
                data_feats_aas_mut.index=data_fit.index
                data_feats_aas_ref.index=data_fit.index
                data_all=pd.concat([data_fit,data_feats_prt,data_feats_aas_mut,data_feats_aas_ref],axis=1)
            else :
                logging.warning("len(data_feats_aas_ref)(%d) != len(data_fit)(%d)" % (len(data_feats_aas_ref),len(data_fit)))    
        else :
            logging.warning("len(data_feats_aas_mut)(%d) != len(data_fit)(%d)" % (len(data_feats_aas_mut),len(data_fit)))    
    else :
        logging.warning("len(data_feats_prt)(%d) != len(data_fit)(%d)" % (len(data_feats_prt),len(data_fit)))    


    data_all["Mutant amino acid"]   =data_all["mut"]
    data_all["Reference amino acid"]=data_all["ref"]
    X_cols=['Mutant amino acid','Reference amino acid'] # X_cols=['mut','ref']
    X_cols=np.array(X_cols+data_feats_prt.columns.tolist() \
                    +data_feats_aas_mut.columns.tolist() \
                    +data_feats_aas_ref.columns.tolist())
    if not keep_mutids:
        data_all=data_all.loc[:,[y_coln]+list(X_cols)]
    elif keep_mutids:
        data_all=data_all.loc[:,["mutids",y_coln]+list(X_cols)]
    return data_all,X_cols

# def y2classes(data_combo,y_coln,mx=None,mn=None,increment=None,quantiles=False):
#     """
#     This converts y values (if numeric) to classes(bins). 
    
#     :param data_combo: dataframe with y and all Xs.
#     :param y_coln: name of column with classes 
#     :param mx: maximum value for binning
#     :param mn: minimum value for binning
#     :param increment: increment for binning
#     :returns data_all: dataframe with (y) classes and all Xs.
#     """
#     if data_combo.applymap(np.isreal).all(0)[y_coln]:
#         if not quantiles:
#             class_bins=zip(np.arange(mn, mx, increment),np.arange(mn+increment, mx+increment, increment))
#         elif quantiles:
#             boundaries_min=[np.percentile(data_combo.loc[~pd.isnull(data_combo.loc[:,y_coln]),y_coln],0),
#             np.percentile(data_combo.loc[~pd.isnull(data_combo.loc[:,y_coln]),y_coln],5),
#             np.percentile(data_combo.loc[~pd.isnull(data_combo.loc[:,y_coln]),y_coln],50)]
#             boundaries_max=boundaries_min[1:]+\
#             [np.percentile(data_combo.loc[~pd.isnull(data_combo.loc[:,y_coln]),y_coln],100)]
#             class_bins=zip(boundaries_min,boundaries_max)
#         for class_bin in class_bins:
#             class_bool=((data_combo.loc[:,y_coln]>class_bin[0]) & (data_combo.loc[:,y_coln]<=class_bin[1]))
#             if sum(class_bool)>100:
#                 data_combo.loc[class_bool, "tmp"]="%.1f < %s <= %.1f" % (class_bin[0],y_coln,class_bin[1])
#             else:
#                 data_combo.loc[class_bool, "tmp"]=np.nan                
# #         data_all=data_combo
#         data_combo.loc[:,y_coln]=data_combo.loc[:,"tmp"]
#     return data_combo

def y2classes(data_combo,y_coln,classes=2):
    if classes==2:
        median=data_combo.loc[:,y_coln].median()
        data_combo.loc[data_combo.loc[:,y_coln]>median,"classes"]="gt median"
        data_combo.loc[data_combo.loc[:,y_coln]<median,"classes"]="lt median"
    return data_combo

def X_cols2numeric(data_all,X_cols,keep_cols=[]):
    """
    This converts features in text form (eg. C, H, L, ..) to numeric (eg. 0, 1, 2, ..)
    
    :param data_all: dataframe with all Xs in columns.
    :param X_cols: name/s of colum/s with feature/s
    :returns data_all: dataframe with all numeric Xs.
    """
    for X_col in X_cols:
        # if not data_all.applymap(np.isreal).all(0)[X_col]:
        if not is_numeric(data_all.loc[:,X_col]):
            if not X_col in keep_cols:
                le = LabelEncoder()
                le.fit(data_all.loc[:,X_col])            
                data_all.loc[:,X_col]=le.transform(data_all.loc[:,X_col])
    return data_all

def X_cols2binary(data,cols=None):    
    if cols==None:
        cols=[]
        for col in data.columns.tolist():
            if not is_numeric(data.loc[:,col]):
                cols.append(col)
    for col in cols:
        classes=list(data.loc[:,col].unique())
        if np.nan in classes:
            classes.remove(np.nan)
        for classi in classes:            
            data.loc[data.loc[:,col]==classi,"%s: %s" % (col,classi)]=1
            data.loc[~(data.loc[:,col]==classi),"%s: %s" % (col,classi)]=0
            data.loc[(data.loc[:,col]==np.nan),"%s: %s" % (col,classi)]=np.nan
        data=data.drop(col,axis=1)
    return data

def binary2classes(y_pred,classes):
    y_pred_classes=[]
    if len(classes)>2:
        for row in y_pred:
            for classi in range(len(classes)):
                if np.sum(row)==1:
                    if row[classi]==1:
                        y_pred_classes.append(classes[classi])
                        break
                else:
                    y_pred_classes.append(np.nan)
                    break
    elif len(classes)==2:
        for row in y_pred:
            for classi in range(len(classes)):
                if row==0:
                    y_pred_classes.append(classes[classi])
                    break
                elif row==1:
                    y_pred_classes.append(classes[classi])
                    break
    return y_pred_classes

def denanrows(data_all):
    """
    This removes rows with any np.nan value/s.  
    
    :param data_all: input dataframe.
    :returns data_all: output dataframe.
    """
    keep_rows_bool=[]
    if "mutids" in data_all.columns.tolist():
        data_all_use=data_all.drop("mutids",axis=1)
    else:
        data_all_use=data_all
    for rowi in range(len(data_all_use)):
        keep_rows_bool.append(all(~pd.isnull(data_all_use.iloc[rowi,:])))
    data_all=data_all.loc[keep_rows_bool,:]
    return data_all

def plot_ROC(y_test,y_score,classes):
    """
    This plots ROC curve.
    
    :param y_test: test split used for predictions.
    :param y_score: probabilities of predictions.
    :param classes: list with unique classes in y
    """
#     fig = plt.figure(figsize=(3,3),dpi=300)#figsize=(11,5))
    fig = plt.figure(figsize=(8.5,3),dpi=300)
    ax_roc = plt.subplot(131)

    mean_tpr = 0.0
    mean_fpr = np.linspace(0, 1, 100)
    classes_to_plot=[]
    if len(classes)>2:
        for classi in range(len(classes)):
#             print np.shape(y_score[0])
#             rows, cols = np.shape(y_score[0])
            fpr, tpr, _ = roc_curve(y_test[:, classi], y_score[0][:, 1])
            if auc(fpr,tpr)<0.5:
                fpr, tpr, _ = roc_curve(y_test[:, classi], y_score[0][:, 0])
#             np.savetxt(classes[classi], np.hstack((fpr,tpr)), delimiter=',')
            if len(np.unique(tpr))>10:
                mean_tpr += np.interp(mean_fpr, fpr, tpr)
                mean_tpr[0] = 0.0
                ax_roc.plot(fpr, tpr, label="%s (AUC=%.2f)" % (classes[classi],auc(fpr,tpr)))
                classes_to_plot.append(classi)
#             else:
#                 logging.warning("error in yscore")
        if len(classes_to_plot)!=0:
            mean_tpr /= len(classes_to_plot)
            mean_tpr[-1] = 1.0
            ax_roc.plot(mean_fpr, mean_tpr,'k', lw=4, label="%s (AUC=%.2f)" % ("mean",auc(mean_fpr,mean_tpr)))
    else:
        fpr, tpr, _ = roc_curve(y_test, y_score[:, 1])
        ax_roc.plot(fpr, tpr, label="%s (AUC=%.2f)" % ("mean",auc(fpr,tpr)))

    ax_roc.plot([0, 1], [0, 1], 'k--')
    ax_roc.set_xlabel("FPR")
    ax_roc.set_ylabel("TPR")
#     ax_roc.legend(loc='lower right',prop={'size':9})
# #     ax_roc.legend(loc='center right', bbox_to_anchor=(2, 0.5))
    ax_roc.legend(loc='lower right', bbox_to_anchor=(2.5, 0.5))
    plt.tight_layout()
    
def plot_importances(importances,X_cols):
    """
    This plots relative importances of features.
    
    :param importances: relative importances of features from `.best_estimator_.feature_importances_`
    :param X_cols: list of features in same order as `importances`
    :returns feature_importances: dataframe with relative importances of features.
    """
    feature_importances=pd.DataFrame(columns=['feature','importances'])
    feature_importances.loc[:,'feature']=X_cols#[indices]
    feature_importances.loc[:,'importances']=importances#[indices]
    fig = plt.figure(figsize=(8,len(feature_importances)*0.25))#figsize=(11,5))
    ax_imp = plt.subplot(1,1,1)

    feature_importances=feature_importances.sort_values("importances",axis=0,ascending=False)
    feature_importances.plot(kind="barh",ax=ax_imp,legend=False)
    ax_imp.invert_yaxis()
    ax_imp.set_yticklabels(feature_importances.loc[:,"feature"])
    ax_imp.set_xlabel("Relative importances")
    ax_imp.set_xticklabels(ax_imp.get_xticks(),rotation=90)
    plt.tight_layout()
    return feature_importances

def run_RF(data_all,X_cols,y_coln,plot_fh="test",test_size=.5):
    """
    This implements Random Forest classifier.
    
    :param data_all: dataframe with columns with features(Xs) and classes(y).
    :param X_cols: list of column names with features.
    :param y_coln: column name of column with classes.
    :param plot_fh: path to output plot file. 
    :returns grid_search: trained classifier object.
    :returns y_test: classes used for testing classifier. 
    :returns y_pred: predicted classes.
    :returns y_score: scores of predicted classes used to plot ROC curve.
    :returns feature_importances: relative importances of features (dataframe).
    """    
    X=data_all.loc[:,list(X_cols)]
    X=X.as_matrix()
    y=data_all.loc[:,y_coln]
    classes=y.unique()
    if len(classes)>1:
        # y.to_csv("testy")
        # print classes
        y = label_binarize(y, classes=classes)

        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=test_size,
                                                            random_state=0)

        model = RandomForestClassifier()
        param_grid = {"n_estimators": [1000],
                      "max_features": ['sqrt', 'auto','log2'],
                      "criterion": ["gini", "entropy"]}

        grid_search = GridSearchCV(model, param_grid=param_grid,cv=5)
        grid_search.fit(X_train,y_train)
        y_pred=grid_search.predict(X_test)
        y_score=grid_search.predict_proba(X_test)
        plot_ROC(y_test,y_score,classes)
        plt.savefig(plot_fh+".pdf",format='pdf')
        plt.savefig(plot_fh);plt.clf();plt.close()
        # model.fit(X_train, y_train)
        # y_score = model.predict_proba(X_test)
        # y_pred = model.predict(X_test)
    #         logging.info("accuracy_score: %.2f, f1_score: %.2f, AUC: %.2f" % (accuracy_score(y_test,y_pred),f1_score(y_test,y_pred),roc_auc_score(y_test,y_pred)))
        importances = grid_search.best_estimator_.feature_importances_
        feature_importances=plot_importances(importances,X_cols)
        plt.savefig(plot_fh.replace("_roc_","_relative_importances_")+".pdf",format='pdf')
        plt.savefig(plot_fh.replace("_roc_","_relative_importances_"));plt.clf();plt.close()
        return grid_search,y_test,y_pred,y_score,feature_importances
    else:
        logging.info("skipping: coz only one class")
                     
def data_fit2ml(data_fit_key,prj_dh,data_feats,y_coln):
    """
    This runs the submodules to run classifier from fitness data (`data_fit`).
    
    :param data_fit_key: in the form <data_fit>/<aas/cds>/<name of file>.
    :param prj_dh: path to project directory.
    :param data_feats: dataframe with features.
    :param y_coln: column name of column with classes (ys). 
    """
    type_form='aas' # in mut_types_form:
    plot_dh="%s/plots/%s" % (prj_dh,type_form)
    data_dh="%s/data_ml/%s" % (prj_dh,type_form)
    plot_fh="%s/fig_ml_roc_%s.png" % (plot_dh,data_fit_key.replace('/','_'))
    data_fh="%s/data_ml_%s" % (data_dh,data_fit_key.replace('/','_'))
    if not exists(plot_fh):
        data_fit=pd.read_csv("%s/%s" % (prj_dh,data_fit_key))
        if np.sum(~data_fit.loc[:,y_coln].isnull())>10:
            if len(data_fit.loc[~data_fit.loc[:,y_coln].isnull(),y_coln].unique())>=3:        
                logging.info("processing: %s" % data_fit_key)
                if not exists(data_fh):
                    if not y_coln=="classes":
                        data_fit=y2classes(data_fit,y_coln)
                        y_coln="classes"
                    data_combo,X_cols=data_fit_feats2combo(data_fit,data_feats,y_coln)
                    # data_all=y2classes(data_combo,y_coln,quantiles=True)# make classes
                    # data_all=y2classes(data_combo,y_coln,\
                    #                    mx=data_combo.loc[:,y_coln].max(),\
                    #                    mn=data_combo.loc[:,y_coln].min(),\
                    #                    increment=data_combo.loc[:,y_coln].median()-data_combo.loc[:,y_coln].min())# make classes
                    data_all=X_cols2numeric(data_combo,X_cols)
                    data_all=denanrows(data_all)# remove nan rows
                    data_all.reset_index().to_csv(data_fh,index=False)
                else:
                    data_all=pd.read_csv(data_fh)
                    X_cols=data_all.columns.tolist()
                    if "class_fit" in X_cols:
                        X_cols.remove("class_fit")
                    if "aasi" in X_cols:
                        X_cols.remove("aasi")
                logging.info("number of samples(mutants) = %d" % len(data_all))
                grid_search,y_test,y_pred,y_score,feature_importances=run_RF(data_all,X_cols,y_coln,plot_fh) #
                try:
                    feature_importances.reset_index().to_csv(data_fh.replace("_ml_","_ml_rel_imp_"),index=False)
                    run_RF(data_all,feature_importances.loc[:,'feature'].head(40).tolist(),\
                           y_coln,plot_fh.replace("_roc_","_roc_top40_"))
                except:
                    logging.info("exception")
            else:
                logging.info("skipping: requires more unique classes")
        else:
            logging.info("skipping %s: requires more samples %d<10" %\
                            (data_fit_key,np.sum(~data_fit.loc[:,y_coln].isnull())))
