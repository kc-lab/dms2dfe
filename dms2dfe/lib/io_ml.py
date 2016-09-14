#!usr/bin/python

# Copyright 2016, Rohan Dandage <rraadd_8@hotmail.com,rohan@igib.in>
# This program is distributed under General Public License v. 3.  

"""
================================
``io_ml``
================================
"""
from os.path import abspath,dirname,exists,basename

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
logging.basicConfig(format='[%(asctime)s] %(levelname)s\tfrom %(filename)s in %(funcName)s(..):%(lineno)d: %(message)s',level=logging.DEBUG) # filename=cfg_xls_fh+'.log'


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

def y2classes(data_combo,y_coln,classes=2):
    if classes==2:
        median=data_combo.loc[:,y_coln].median()
        data_combo.loc[data_combo.loc[:,y_coln]>=median,"classes"]="gt median"
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

def zscore(df,col):
    df.loc[:,col] = (df.loc[:,col]-df.loc[~pd.isnull(df.loc[:,col]),col].mean())/df.loc[~pd.isnull(df.loc[:,col]),col].std()
    return df

def rescalecols(data_combo,kind="zscore"):
    for col in data_combo.columns.tolist():
#         print col
        if is_numeric(data_combo.loc[:,col]):
            if kind=="zscore" and data_combo.loc[~pd.isnull(data_combo.loc[:,col]),col].std() !=0:
#                 print col
                data_combo=zscore(data_combo,col)
            else:
                data_combo.loc[:,col]=data_combo.loc[:,col]
    return data_combo

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

def denanrows(data_all,condi="any"):
    """
    This removes rows with any np.nan value/s.  
    
    :param data_all: input dataframe.
    :param condi: conditions for deletion of rows ["any": if any element is nan ,default | "all" : if all elements are nan] 
    :returns data_all: output dataframe.
    """
    keep_rows_bool=[]
    if "mutids" in data_all.columns.tolist():
        data_all_use=data_all.drop("mutids",axis=1)
    else:
        data_all_use=data_all
    for rowi in range(len(data_all_use)):
        if condi=="any":
            keep_rows_bool.append(all(~pd.isnull(data_all_use.iloc[rowi,:])))
        if condi=="all":
            keep_rows_bool.append(any(~pd.isnull(data_all_use.iloc[rowi,:])))
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
            logging.info("mean AUC = %.2f" % auc(fpr,tpr))
    else:
        fpr, tpr, _ = roc_curve(y_test, y_score[:, 1])
        ax_roc.plot(fpr, tpr, label="%s (AUC=%.2f)" % ("mean",auc(fpr,tpr)))
        logging.info("mean AUC = %.2f" % auc(fpr,tpr))

    ax_roc.plot([0, 1], [0, 1], 'k--')
    ax_roc.set_xlabel("FPR")
    ax_roc.set_ylabel("TPR")
#     ax_roc.legend(loc='lower right',prop={'size':9})
# #     ax_roc.legend(loc='center right', bbox_to_anchor=(2, 0.5))
    ax_roc.legend(loc='lower right',
                  # bbox_to_anchor=(2.5, 0.5)
                 )
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

def run_RF(data_all,X_cols,y_coln,plot_fh="test",test_size=.5,data_test=None):
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
    from sklearn.ensemble import RandomForestClassifier

    X=data_all.loc[:,list(X_cols)]
    X=X.as_matrix()

    y=data_all.loc[:,y_coln]
    classes=y.unique()
    y=y.as_matrix()
    y = label_binarize(y, classes=classes)
    if len(classes)==2:
        y=np.array([i[0] for i in y])

    if len(classes)>1:
        if test_size!=0:
            X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=test_size,
                                                            random_state=0)
        else :
            X_train=X
            y_train=y
            X_test_df=data_test.loc[:,list(X_cols)]
#             X_test_df.to_csv("test_xtest")
            X_test_df=denanrows(X_test_df)
            X_test=X_test_df.as_matrix()
            y_test=None

        model = RandomForestClassifier(random_state =88)
        param_grid = {"n_estimators": [1000],
                      "max_features": [None,'sqrt','log2'],
                      "min_samples_leaf":[1,25,50,100],
                      "criterion": ["gini", "entropy"]}

        grid_search = GridSearchCV(model, param_grid=param_grid,cv=5)
#     return grid_search,X_train,y_train
#         np.savetxt("X_train", X_train,delimiter=",")
#         np.savetxt("y_train", y_train,delimiter=",")
        grid_search.fit(X_train,y_train)

        y_pred=grid_search.predict(X_test)
        y_score=grid_search.predict_proba(X_test)

        if test_size!=0:    
            data_preds=None
            plot_ROC(y_test,y_score,classes)
            plt.savefig(plot_fh+".pdf",format='pdf')
            # plt.savefig(plot_fh);plt.clf();plt.close()
        else:
            data_preds=X_test_df
#             np.savetxt("y_pred", y_pred,delimiter=",")
#             print classes
            data_preds[y_coln]=binary2classes(y_pred,classes)

        importances = grid_search.best_estimator_.feature_importances_
        feature_importances=plot_importances(importances,X_cols)
        plt.savefig(plot_fh.replace("_roc_","_relative_importances_")+".pdf",format='pdf')
        # plt.savefig(plot_fh.replace("_roc_","_relative_importances_"));plt.clf();plt.close()
    return grid_search,y_test,y_pred,y_score,feature_importances,data_preds

def run_RF_regress(data_all,X_cols,y_coln,test_size=0.5,data_test=None,plot_fh="test"):
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
    from sklearn.ensemble import RandomForestRegressor

    X=data_all.loc[:,list(X_cols)]
    X=X.as_matrix()
    y=data_all.loc[:,y_coln]
    y=y.as_matrix()
    
    if test_size!=0:
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=test_size,
                                                        random_state=0)
    else :
        X_train=X
        y_train=y
        X_test=data_test.loc[:,list(X_cols)].as_matrix()
        y_test=None

    model = RandomForestRegressor(random_state =88)
    param_grid = {"n_estimators": [1000],
                  "max_features": [None,'sqrt','log2'],
                  "min_samples_leaf":[1,25,50,100],
                  "criterion": ["mse"],
                  "oob_score": [True],
                 }

    grid_search = GridSearchCV(model, param_grid=param_grid,cv=5)
    grid_search.fit(X_train,y_train)
    y_pred=grid_search.predict(X_test)

    if test_size!=0:    
        data_preds=None
        # print grid_search.score(X_test, y_test)
    else:
        data_preds=data_test.loc[:,list(X_cols)]
        data_preds[y_coln]=y_pred


    importances = grid_search.best_estimator_.feature_importances_
    feature_importances=plot_importances(importances,X_cols)
    plt.savefig(plot_fh.replace("_roc_","_relative_importances_")+".pdf",format='pdf')
    # plt.savefig(plot_fh.replace("_roc_","_relative_importances_"));plt.clf();plt.close()
    return grid_search,y_test,y_pred,feature_importances,data_preds

def data_fit2ml(data_fit_key,prj_dh,data_feats):
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
    y_coln_classi="FCA"
    if not exists(data_fh.replace("data_ml_","data_ml_regress_preds_")):
        data_fit=pd.read_csv("%s/%s" % (prj_dh,data_fit_key))
        if np.sum(~data_fit.loc[:,y_coln_classi].isnull())>10:
            # if len(data_fit.loc[~data_fit.loc[:,y_coln_classi].isnull(),y_coln_classi].unique())>=3:        
            logging.info("processing: %s" % data_fit_key)
            if not exists(data_fh.replace("data_ml_","data_ml_classi_feature_importances_")):
                if not exists(data_fh.replace("data_ml_","data_ml_classi_train_")):
                    data_fit=y2classes(data_fit,y_coln_classi)
                    y_coln_classi="classes"
                    data_combo,X_cols_classi=data_fit_feats2combo(data_fit,data_feats,y_coln_classi,keep_mutids=True)
                    data_combo=data_combo.set_index("mutids",drop=True)
                    data_ml=X_cols2binary(data_combo.drop(y_coln_classi,axis=1))
                    data_ml.loc[:,y_coln_classi]=data_combo.loc[:,y_coln_classi]
                    data_ml=rescalecols(data_ml)
                    data_classi_train=denanrows(data_ml)                    
#                         data_classi_train=y2classes(data_classi_train,y_coln_classi)
                    data_ml_mutids=list(data_ml.index.values)                
                    data_classi_train_mutids=list(data_classi_train.index.values)
                    data_classi_test_mutids=[mutid for mutid in data_ml_mutids if not mutid in data_classi_train_mutids]
                    data_classi_test=data_ml.loc[data_classi_test_mutids,:]

                    data_combo.reset_index().to_csv(data_fh.replace("data_ml_","data_ml_combo_"),index=False)
                    data_ml.reset_index().to_csv(data_fh.replace("data_ml_","data_ml_classi_all_"),index=False)
                    data_classi_train.reset_index().to_csv(data_fh.replace("data_ml_","data_ml_classi_train_"),index=False)
                    data_classi_test.reset_index().to_csv(data_fh.replace("data_ml_","data_ml_classi_test_"),index=False)
                else:                
                    data_classi_train=pd.read_csv(data_fh.replace("data_ml_","data_ml_classi_train_"))
                    data_classi_test=pd.read_csv(data_fh.replace("data_ml_","data_ml_classi_test_"))

                    data_classi_train  =data_classi_train.set_index("mutids",drop=True)
                    data_classi_test  =data_classi_test.set_index("mutids",drop=True)
                    y_coln_classi="classes"
                logging.info("number of mutants used for training = %d" % len(data_classi_train))
                logging.info("this step would take a 10 to 15min to complete.")

                X_cols_classi=data_classi_train.columns.tolist()
                X_cols_classi.remove(y_coln_classi)
#                     # 0.5 split for check up
#                     grid_search_classi,y_test_classi,y_pred_classi,y_score_classi,feature_importances_classi1,data_preds_classi\
#                     =run_RF(data_classi_train,X_cols_classi,y_coln_classi,plot_fh=plot_fh.replace("fig_ml_","fig_ml_classi1_"),\
#                         test_size=0.5,data_test=data_classi_test) #
#                     # no split 
#                     grid_search_classi,y_test_classi,y_pred_classi,y_score_classi,feature_importances_classi2,data_preds_classi\
#                     =run_RF(data_classi_train,feature_importances_classi1.loc[feature_importances_classi1.loc[:,"importances"]>0.002,'feature'].tolist(),\
#                             y_coln_classi,plot_fh=plot_fh.replace("fig_ml_","fig_ml_classi2_"),\
#                             test_size=0,data_test=data_classi_test) #

#                     feature_importances_classi1.to_csv(data_fh.replace("data_ml_","data_ml_classi1_feature_importances_"))
#                     feature_importances_classi2.to_csv(data_fh.replace("data_ml_","data_ml_classi2_feature_importances_"))
                # classi
                grid_search_classi,y_test_classi,y_pred_classi,y_score_classi,feature_importances_classi,data_preds_classi\
                =run_RF(data_classi_train,X_cols_classi,y_coln_classi,plot_fh=plot_fh.replace("fig_ml_","fig_ml_classi1_"),\
                    test_size=0.5,data_test=data_classi_test) #

                feature_importances_classi.to_csv(data_fh.replace("data_ml_","data_ml_classi_feature_importances_"))
#                     data_preds_classi.to_csv(data_fh.replace("data_ml_","data_ml_classi_preds_"))
#                     data_classi_all=data_classi_train.append(data_preds_classi)
#                     data_classi_all.to_csv(data_fh.replace("data_ml_","data_ml_classi_all_"))
            else:
                y_coln_classi="classes"
#                     feature_importances_classi2=pd.read_csv(data_fh.replace("data_ml_","data_ml_classi_feature_importances_"))
                feature_importances_classi=pd.read_csv(data_fh.replace("data_ml_","data_ml_classi_feature_importances_"))
#                     data_regress_test=pd.read_csv(data_fh.replace("data_ml_","data_ml_classi_preds_"))
                data_classi_test=pd.read_csv(data_fh.replace("data_ml_","data_ml_classi_test_"))
                data_classi_test  =data_classi_test.set_index("mutids",drop=True)
                data_classi_train=pd.read_csv(data_fh.replace("data_ml_","data_ml_classi_train_"))
                data_classi_train  =data_classi_train.set_index("mutids",drop=True)
            #regress
            y_coln_regress="FCA"
            data_regress_train=data_classi_train                
            data_regress_train=X_cols2binary(data_regress_train,[y_coln_classi])
            data_regress_train.loc[:,y_coln_regress]\
            =data_fit.set_index("mutids").loc[data_classi_train.index.values,y_coln_regress]
            data_regress_train=denanrows(data_regress_train)
            data_regress_test=data_classi_test                
            if y_coln_classi in data_regress_test.columns.tolist(): 
                data_regress_test=data_regress_test.drop(y_coln_classi,axis=1)
            else:
                data_regress_test=denanrows(data_regress_test)
                data_regress_test=X_cols2binary(data_regress_test,[y_coln_classi])
            data_regress_test=denanrows(data_regress_test)

#                 X_cols_regress_class_fit=[col for col in data_regress_train.columns.tolist() if y_coln_classi in col]
            X_cols_regress_class_fit=[]
#                 X_cols_regress=[col for col in X_cols_regress_class_fit+list(feature_importances_classi2.loc[:,'feature']) \
#                                 if col in data_regress_test.columns.tolist()]
            X_cols_regress=[col for col in X_cols_regress_class_fit+list(feature_importances_classi.loc[feature_importances_classi.loc[:,"importances"]>0.002,'feature']) \
                            if col in data_regress_test.columns.tolist()]
            if y_coln_regress in X_cols_regress:
                X_cols_regress.remove(y_coln_regress)            
            data_regress_train.to_csv(data_fh.replace("data_ml_","data_ml_regress_train_"))
            data_regress_test.to_csv(data_fh.replace("data_ml_","data_ml_regress_test_"))
            try:
                grid_search_regress,y_test_regress,y_pred_regress,feature_importances_regress,data_preds_regress\
                =run_RF_regress(data_regress_train,X_cols_regress,y_coln_regress,\
    test_size=0,data_test=data_regress_test,plot_fh=plot_fh.replace("fig_ml_","fig_ml_regress_"))

                data_preds_regress.to_csv(data_fh.replace("data_ml_","data_ml_regress_preds_"))
                data_regress_all=data_preds_regress.append(data_regress_train)
                data_regress_all.to_csv(data_fh.replace("data_ml_","data_ml_regress_all_"))
                data_regress2data_fit(prj_dh,data_fit_key,data_regress_all)
                return data_regress_all
            except:
                logging.info("skipping: %s : requires more data" % basename(data_fit_key))
        else:
            logging.info("skipping %s: requires more samples %d<10" %\
                            (data_fit_key,np.sum(~data_fit.loc[:,y_coln].isnull())))
    else:
        data_regress_all=pd.read_csv(data_fh.replace("data_ml_","data_ml_regress_all_"))
        data_regress2data_fit(prj_dh,data_fit_key,data_regress_all)
        return data_regress_all
    
def data_regress2data_fit(prj_dh,data_fit_key,data_regress_all):
    from dms2dfe.lib.io_nums import str2num
    from dms2dfe.lib.io_mut_files import rescale_fitnessbysynonymous,class_fit

    data_fit=pd.read_csv("%s/%s" % (prj_dh,data_fit_key))
    data_fit=data_fit.loc[:,["mutids","FCA"]].set_index("mutids",drop=True)
    data_fit_infered=data_regress_all.reset_index().loc[:,["mutids","FCA"]].set_index("mutids",drop=True)
    for mutid in data_fit.index.values:
        if not mutid in data_fit_infered.index.values:
            data_fit_infered.loc[mutid,"FCA"]=data_fit.loc[mutid,"FCA"]
    
    data_fit_infered=pd.DataFrame(data_fit_infered.reset_index())
    # print data_fit_infered.columns
    # str2num()
    data_fit_infered.loc[:,'refi']=[str2num(mutid) for mutid in data_fit_infered.loc[:,"mutids"].tolist()]
    data_fit_infered.loc[:,'ref']=[mutid[0] for mutid in data_fit_infered.loc[:,"mutids"].tolist()]
    data_fit_infered.loc[:,'mut']=[mutid[-1] for mutid in data_fit_infered.loc[:,"mutids"].tolist()]
    data_fit_infered.loc[:,'refrefi']=[("%s%03d" % (mutid[0],str2num(mutid))) for mutid in data_fit_infered.loc[:,"mutids"].tolist()]
    data_fit_infered.loc[:,'FCS']=data_fit_infered.loc[(data_fit_infered.loc[:,'ref']==data_fit_infered.loc[:,'mut']),'FCA']
    # data_fit_infered.head()
    # data_fit_infered.to_csv("data_fit_infered")
    data_fit_infered=rescale_fitnessbysynonymous(data_fit_infered)
    data_fit_infered=class_fit(data_fit_infered)
    data_fit_infered.loc[:,'FiS']=\
    data_fit_infered.loc[(data_fit_infered.loc[:,'ref']==data_fit_infered.loc[:,'mut']),'FiA']
    data_fit_infered=data_fit_infered.sort_values(by="refi",axis=0)
    data_fit_infered.to_csv("%s/%s_inferred" % (prj_dh,data_fit_key))
    
    
    