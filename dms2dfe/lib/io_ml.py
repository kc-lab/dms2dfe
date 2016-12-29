#!usr/bin/python

# Copyright 2016, Rohan Dandage <rraadd_8@hotmail.com,rohan@igib.in>
# This program is distributed under General Public License v. 3.  

"""
================================
``io_ml``
================================
"""
from os.path import abspath,dirname,exists,basename

from sklearn.cross_validation import train_test_split,KFold
from sklearn.preprocessing import LabelEncoder,label_binarize
from sklearn.metrics import roc_curve, auc
from sklearn.grid_search import GridSearchCV
from sklearn.metrics import confusion_matrix,classification_report,regression

from dms2dfe.lib.io_data_files import read_pkl,to_pkl
from dms2dfe.lib.io_dfs import set_index
from dms2dfe.lib.io_nums import is_numeric
from dms2dfe.lib.io_plots import saveplot,get_axlims

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg') # no Xwindows

import matplotlib.pyplot as plt
import seaborn as sns

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

def y2classes(data_combo,y_coln,classes=2,
             middle_percentile_skipped=0):
    data_combo.loc[:,'classes']=np.nan
    if classes==2:
        median=data_combo.loc[:,y_coln].median()
        if middle_percentile_skipped==0:            
            data_combo.loc[data_combo.loc[:,y_coln]>=median,"classes"]="high"
            data_combo.loc[data_combo.loc[:,y_coln]<median,"classes"]="low"
        else:
            up_bound=data_combo.loc[~pd.isnull(data_combo.loc[:,y_coln]),y_coln].quantile(0.5+middle_percentile_skipped/2)            
            lw_bound=data_combo.loc[~pd.isnull(data_combo.loc[:,y_coln]),y_coln].quantile(0.5-middle_percentile_skipped/2)
            # print up_bound
            # print lw_bound            
            data_combo.loc[data_combo.loc[:,y_coln]>up_bound,"classes"]="high"
            data_combo.loc[data_combo.loc[:,y_coln]<lw_bound,"classes"]="low"
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
        data_all_use=data_all.copy()
    for rowi in range(len(data_all_use)):
        if condi=="any":
            keep_rows_bool.append(all(~pd.isnull(data_all_use.iloc[rowi,:])))
        if condi=="all":
            keep_rows_bool.append(any(~pd.isnull(data_all_use.iloc[rowi,:])))
    data_all=data_all.loc[keep_rows_bool,:]
    return data_all

def denan(data_all,axis,condi="any"):
    """
    This removes rows with any np.nan value/s.  
    
    :param data_all: input dataframe.
    :param condi: conditions for deletion of rows ["any": if any element is nan ,default | "all" : if all elements are nan] 
    :returns data_all: output dataframe.
    """
    import logging
    logging.info("denan: original: rows=%s cols=%s" % data_all.shape)

    if axis=='both':
        condi_cols=condi.split(' ')[0]
        condi_rows=condi.split(' ')[1]
    if axis=='rows' or axis==0:
        condi_rows=condi        
    if axis=='cols' or axis==1:
        condi_cols=condi
    if axis=='cols' or axis==1 or axis=='both':
        data_all_use=data_all.copy()
        keep_bool=[]
        for col in data_all_use:
            if condi_cols=="any":
                keep_bool.append(all(~pd.isnull(data_all_use.loc[:,col])))
            if condi_cols=="all":
                keep_bool.append(any(~pd.isnull(data_all_use.loc[:,col])))
        data_all=data_all.loc[:,keep_bool]
        logging.info("denan: cols:      rows=%s cols=%s" % data_all.shape)
    if axis=='rows' or axis==0 or axis=='both':
        data_all_use=data_all.copy()
        keep_bool=[]
        for rowi in range(len(data_all_use)):
            if condi_rows=="any":
                keep_bool.append(all(~pd.isnull(data_all_use.iloc[rowi,:])))
            if condi_rows=="all":
                keep_bool.append(any(~pd.isnull(data_all_use.iloc[rowi,:])))
        data_all=data_all.loc[keep_bool,:]
        logging.info("denan: rows:      rows=%s cols=%s" % data_all.shape)
    return data_all

def plot_ROC(y_test,y_score,classes,lw=2,
             ax_roc=None,annotate=True,
            get_auc=False,reference_line=True,plot_fh=None):
    """
    This plots ROC curve.
    
    :param y_test: test split used for predictions.
    :param y_score: probabilities of predictions.
    :param classes: list with unique classes in y
    """
    # fig = 
    if ax_roc is None: 
        plt.figure(figsize=(3,3),dpi=300)#figsize=(11,5))
        ax_roc = plt.subplot(111)
        
    mean_tpr = 0.0
    mean_fpr = np.linspace(0, 1, 100)
    classes_to_plot=[]
    if len(classes)>2:
        for classi in range(len(classes)):
            fpr, tpr, _ = roc_curve(y_test[:, classi], y_score[0][:, 1])
            if auc(fpr,tpr)<0.5:
                fpr, tpr, _ = roc_curve(y_test[:, classi], y_score[0][:, 0])
            if len(np.unique(tpr))>10:
                mean_tpr += np.interp(mean_fpr, fpr, tpr)
                mean_tpr[0] = 0.0
                ax_roc.plot(fpr, tpr,lw=lw)#, label="%s (AUC=%.2f)" % (classes[classi],auc(fpr,tpr)))
                ax_roc.annotate("%s (AUC=%.2f)" % (classes[classi],auc(fpr,tpr)), xy=(2, 1), xytext=(2, 1))
                classes_to_plot.append(classi)
        if len(classes_to_plot)!=0:
            mean_tpr /= len(classes_to_plot)
            mean_tpr[-1] = 1.0
            ax_roc.plot(mean_fpr, mean_tpr, lw=lw)#, label="%s (AUC=%.2f)" % ("mean",auc(mean_fpr,mean_tpr)))
            logging.info("mean AUC = %.2f" % auc(fpr,tpr))
    else:
        fpr, tpr, _ = roc_curve(y_test, y_score[:, 1])
        if auc(fpr,tpr)<0.5:
            fpr, tpr, _ = roc_curve(y_test, y_score[:, 0])
            logging.info("mean AUC = %.2f" % (1-auc(fpr,tpr)))
        else:
            logging.info("mean AUC = %.2f" % auc(fpr,tpr))
        ax_roc.plot(fpr, tpr,lw=lw)#, label="%s (AUC=%.2f)" % ("mean",auc(fpr,tpr)))

    # get auc score
    if auc(fpr,tpr)<0.5:
        auc_score=1-auc(fpr,tpr)
    else:
        auc_score=auc(fpr,tpr)

    if annotate:
        ax_roc.annotate("AUC = %.2f" % auc_score, 
                    xy=(0.45, 0), xytext=(0.45, 0))
    if reference_line:
        ax_roc.plot([0, 1], [0, 1], 'k--')
    ax_roc.set_xlim([-0.01,1.01])
    ax_roc.set_ylim([-0.01,1.01])
    ax_roc.set_xlabel("FPR")
    ax_roc.set_ylabel("TPR")
#     ax_roc.legend(loc='lower right',
#                  )
#     ax_roc.grid(color='lightgrey')
#     plt.axis('equal')
    # plt.tight_layout()
    saveplot(plot_fh,transparent=False)
    if get_auc:
        return auc_score
    else:
        return ax_roc       
def plot_importances(feature_importances,plot_fh=None,data_out_fh=None):
    """
    This plots relative importances of features.
    
    :param importances: relative importances of features from `.best_estimator_.feature_importances_`
    :param X_cols: list of features in same order as `importances`
    :returns feature_importances: dataframe with relative importances of features.
    """
    fig = plt.figure(figsize=(8,len(feature_importances)*0.25))#figsize=(11,5))
    ax_imp = plt.subplot(1,1,1)

    feature_importances=feature_importances.sort_values("Importance",axis=0,ascending=False)
    feature_importances.plot(kind="barh",ax=ax_imp,legend=False)
    ax_imp.invert_yaxis()
    ax_imp.set_yticklabels(feature_importances.loc[:,"Feature"])
    ax_imp.set_xlabel("Feature importance")
    ax_imp.set_xticklabels(ax_imp.get_xticks(),rotation=90)
    ax_imp.grid()
    saveplot(plot_fh)
    if not data_out_fh is None:
        feature_importances.to_csv(data_out_fh)
    return feature_importances

def get_RF_ci(RF_type,RF_classi,X_train,X_test,y_test,y_score,
                classes=['yes','no'],plot_fh=None):
    import forestci as fci
    # calculate inbag and unbiased variance
    inbag = fci.calc_inbag(X_train.shape[0], RF_classi)
    V_IJ_unbiased = fci.random_forest_error(RF_classi,inbag, X_train,
                                                 X_test)
    # Plot forest prediction for emails and standard deviation for estimates
    # Blue points are spam emails; Green points are non-spam emails
    idx = np.where(y_test == 1)[0]
    fig=plt.figure(figsize=[3,3])
    ax=plt.subplot(111)
    if RF_type=='classi':
        ax.errorbar(y_score[idx, 1], np.sqrt(V_IJ_unbiased[idx]),
                     fmt='.', alpha=0.75, label=classes[0])

        idx = np.where(y_test == 0)[0]
        ax.errorbar(y_score[idx, 1], np.sqrt(V_IJ_unbiased[idx]),
                     fmt='.', alpha=0.75, label=classes[1])

        ax.set_xlabel('Prediction probability')
        ax.set_ylabel('Standard deviation')
        space=0.3
        ax.set_ylim([ax.get_ylim()[0]*(1+space),
                     ax.get_ylim()[1]*(1+space)])
        leg=ax.legend(loc='upper right',frameon=True)
        leg.get_frame().set_alpha(0.5)
        # plt.axis('equal')
    if RF_type=='regress':
        # Plot error bars for predicted MPG using unbiased variance
        ax.errorbar(y_test, y_score, yerr=np.sqrt(V_IJ_unbiased), fmt='o')
        xlim,ylim=get_axlims(y_test,y_score,
                             space=0.1,equal=True)
        ax.plot(xlim,xlim, '--',color='gray')
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_xlabel('Test')
        ax.set_ylabel('Predicted')
        results="$R^{2}$=%0.2f\nRMSE=%0.2f" % (np.sqrt(regression.r2_score(y_test,y_score)),
             regression.mean_absolute_error(y_test,y_score))
        ax.text(0, 1, results,
            horizontalalignment='left',
            verticalalignment='top',
            transform=ax.transAxes)
    ax.grid(True)
    saveplot(plot_fh)
    
def get_RF_cm(y_test, y_pred,classes,plot_fh=None,data_out_fh=None):
    fig=plt.figure(figsize=[2.5,2])
    data=pd.DataFrame(confusion_matrix(y_test, y_pred))
    data.columns=classes
    data.index=classes
    ax=sns.heatmap(data,cmap='bone_r')
#     plt.axis('equal')
    saveplot(plot_fh)
    if not data_out_fh is None:
        data.to_csv(data_out_fh)

def get_RF_cr(y_test,y_pred,classes,data_out_fh=None):
    s=classification_report(y_test, y_pred)
    for i,line in enumerate(s.split('\n')):
        line=line.replace(' / ','/')
        if not line=='':
            if i==0:
                cols=line.split()
                data=pd.DataFrame(columns=cols)
            else:
                for coli,col in enumerate(cols):
                    data.loc[line.split()[0],cols[coli]]=line.split()[coli+1]            

    data.index=list(classes)+[data.index.tolist()[2]]
    if not data_out_fh is None:
        data.to_csv(data_out_fh)
    return data

def run_RF_classi(data_all,X_cols,y_coln,
           test_size=0.34,data_test=None,data_out_fh=None):
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
                                                            random_state=88)
        else :
            X_train=X
            y_train=y
            X_test_df=data_test.loc[:,list(X_cols)]
            X_test_df=denan(X_test_df,axis='both',condi='all any')
            X_test=X_test_df.as_matrix()
            y_test=None

        model = RandomForestClassifier(random_state =88)
        param_grid = {"n_estimators": [1000],
                      "max_features": ['sqrt'],#[None,'sqrt','log2'],
                      "min_samples_leaf":[1],#[1,25,50,100],
                      "criterion": ['entropy'],#["gini", "entropy"]
                     }

        grid_search = GridSearchCV(model, param_grid=param_grid,cv=10)
        grid_search.fit(X_train,y_train)

        y_pred=grid_search.predict(X_test)
        if test_size!=0:    
            data_preds=None
        else:
            data_preds=X_test_df
            data_preds[y_coln]=binary2classes(y_pred,classes)

    featimps=pd.DataFrame(columns=['Feature','Importance'])
    featimps.loc[:,'Feature']=X_cols#[indices]
    featimps.loc[:,'Importance']=grid_search.best_estimator_.feature_importances_

    data={'RF_classi':grid_search,
          'X_train':X_train,
          'X_test':X_test,
          'y_test':y_test,
          'y_score':grid_search.predict_proba(X_test),
          'classes':classes,
          'X_cols':X_cols,
          'features':X_cols,
          'featimps':featimps,
          'y_pred':y_pred,
         'data_preds':data_preds}
    to_pkl(data,data_out_fh)            
    return grid_search,data_preds
    
def gcv2rfc(gridcv):
    classi=gridcv.best_estimator_
    classi.n_estimators=gridcv.param_grid['n_estimators'][0]
    classi.estimators_=gridcv.best_estimator_.estimators_
    return classi
def get_RF_classi_metrics(data_classi_fh,data_dh='data_ml/',plot_dh='plots/'):
    data_classi=read_pkl(data_classi_fh)
    RF_classi=data_classi['RF_classi']
    X_train=data_classi['X_train']
    X_test=data_classi['X_test']
    y_test=data_classi['y_test']
    y_pred=data_classi['y_pred']
    y_score=data_classi['y_score']
    classes=data_classi['classes']
    featimps=data_classi['featimps']
    #roc
    plot_type='roc'
    plot_fh="%s_%s_.pdf" % (data_classi_fh.replace(data_dh,plot_dh),plot_type)
    plot_ROC(y_test,y_score,classes,plot_fh=plot_fh)
    #featimps
    plot_type='featimps'
    plot_fh="%s_%s_.pdf" % (data_classi_fh.replace(data_dh,plot_dh),plot_type)
    data_out_fh="%s_%s_.csv" % (data_classi_fh,plot_type)
    plot_importances(featimps,plot_fh=plot_fh,data_out_fh=data_out_fh)

    # ci :confidence intervals
    plot_type='ci'
    plot_fh="%s_%s_.pdf" % (data_classi_fh.replace(data_dh,plot_dh),plot_type)
    get_RF_ci('classi',gcv2rfc(RF_classi),X_train,X_test,
                 y_test,y_score,classes=classes,plot_fh=plot_fh)    

    # cm : confusion matrix
    plot_type='cm'
    plot_fh="%s_%s_.pdf" % (data_classi_fh.replace(data_dh,plot_dh),plot_type)
    data_out_fh="%s_%s_.csv" % (data_classi_fh,plot_type)
    get_RF_cm(y_test, y_pred,classes,plot_fh=plot_fh,data_out_fh=data_out_fh)
    #cr : classi report
    plot_type='cr'
    plot_fh="%s_%s_.pdf" % (data_classi_fh.replace(data_dh,plot_dh),plot_type)
    data_out_fh="%s_%s_.csv" % (data_classi_fh,plot_type)
    get_RF_cr(y_test,y_pred,classes,data_out_fh=data_out_fh)

def run_RF_regress(data_all,X_cols,y_coln,
                   test_size=0.5,data_test=None,data_out_fh=None):
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
                                                        random_state=88)
    else :
        X_train=X
        y_train=y
        X_test=data_test.loc[:,list(X_cols)].as_matrix()
        y_test=None

    model = RandomForestRegressor(random_state =88)
    param_grid = {"n_estimators": [3000],#[1000,2000,4000],#
                  "max_features": ['sqrt'],#[None,'sqrt','log2'],
                  "min_samples_leaf":  [1],#[1,25,50,100],
                  "criterion": ["mse"],
                  "oob_score": [True],
                 }

    grid_search = GridSearchCV(model, param_grid=param_grid,cv=10)
    grid_search.fit(X_train,y_train)
    y_pred=grid_search.predict(X_test)

    if test_size!=0:    
        data_preds=None
        # print grid_search.score(X_test, y_test)
    else:
        data_preds=data_test.loc[:,list(X_cols)]
        data_preds[y_coln]=y_pred

    featimps=pd.DataFrame(columns=['Feature','Importance'])
    featimps.loc[:,'Feature']=X_cols#[indices]
    featimps.loc[:,'Importance']=grid_search.best_estimator_.feature_importances_

    data={'RF_regress':grid_search,
          'X_train':X_train,
          'X_test':X_test,
          'y_test':y_test,
          'X_cols':X_cols,
          'features':X_cols,
          'featimps':featimps,
          'y_pred':y_pred,
         'data_preds':data_preds}
    to_pkl(data,data_out_fh)            
    return grid_search,data_preds

def get_RF_regress_metrics(data_regress_fh,data_dh='data_ml/',plot_dh='plots/'):
    data_regress=read_pkl(data_regress_fh)
    RF_regress=data_regress['RF_regress']
    X_train=data_regress['X_train']
    X_test=data_regress['X_test']
    y_test=data_regress['y_test']
    y_pred=data_regress['y_pred']
    featimps=data_regress['featimps']
    #featimps
    plot_type='featimps'
    plot_fh="%s_%s_.pdf" % (data_regress_fh.replace(data_dh,plot_dh),plot_type)
    data_out_fh="%s_%s_.csv" % (data_regress_fh,plot_type)
    importances = RF_regress.best_estimator_.feature_importances_
    feature_importances=plot_importances(featimps,plot_fh=plot_fh,data_out_fh=data_out_fh)

    # ci :confidence intervals
    plot_type='ci'
    plot_fh="%s_%s_.pdf" % (data_regress_fh.replace(data_dh,plot_dh),plot_type)
    get_RF_ci('regress',gcv2rfc(RF_regress),X_train,X_test,
                 y_test,y_pred,plot_fh=plot_fh)    

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
    plot_classi_fh="%s/fig_ml_classi_%s.pdf" % (plot_dh,data_fit_key.replace('/','_'))
    plot_regress_fh="%s/fig_ml_regress_%s.pdf" % (plot_dh,data_fit_key.replace('/','_'))
    data_fh="%s/data_ml_%s" % (data_dh,data_fit_key.replace('/','_'))
    grid_search_classi_fh=data_fh.replace("data_ml_","data_ml_classi_")+'.pkl'
    grid_search_regress_fh=data_fh.replace("data_ml_","data_ml_regress_")+'.pkl'
    grid_search_classi_metrics_fh=data_fh.replace("data_ml_","data_ml_classi_metrics_")+'.pkl'
    grid_search_regress_metrics_fh=data_fh.replace("data_ml_","data_ml_regress_metrics_")+'.pkl'

    y_coln_classi="FCA_norm"
    if not exists(grid_search_regress_fh):
        data_fit_fh="%s/%s" % (prj_dh,data_fit_key)
        data_fit=pd.read_csv(data_fit_fh)
        if np.sum(~data_fit.loc[:,y_coln_classi].isnull())>10:
            logging.info("processing: %s" % data_fit_key)
            if not exists(grid_search_classi_fh):
                if not exists(data_fh.replace("data_ml_","data_ml_classi_train_")):
                    data_fit=y2classes(data_fit,y_coln_classi,
                                       middle_percentile_skipped=0.20)
                    data_ml_mutids=list(data_fit.loc[:,'mutids'])                
                    y_coln_classi="classes"
                    data_fit=set_index(data_fit,"mutids")
                    data_feats=set_index(data_feats,"mutids")
                    X_cols_classi=data_feats.columns.tolist()
                    data_combo=pd.concat([data_feats,
                                          data_fit.loc[:,y_coln_classi]],axis=1)
                    data_combo.index.name='mutids'
                    data_ml=X_cols2binary(data_combo.drop(y_coln_classi,axis=1))
                    data_ml.loc[:,y_coln_classi]=data_combo.loc[:,y_coln_classi]
                    data_ml=rescalecols(data_ml)
                    data_classi_train=denan(data_ml,axis='both',condi='all any')

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
                # classi
                grid_search_classi,data_preds=run_RF_classi(data_classi_train,X_cols_classi,y_coln_classi,
                        test_size=0.34,data_test=data_classi_test,data_out_fh=grid_search_classi_fh) #                     
                get_RF_classi_metrics(grid_search_classi_fh,data_dh='data_ml/',plot_dh='plots/')
                
                # classi metrics
                feature_importances_classi_fh="%s_%s_.csv" % (grid_search_classi_fh,'featimps')
                feature_importances_classi=pd.read_csv(feature_importances_classi_fh)
                X_cols_classi_metrics_selected=feature_importances_classi\
                .sort_values(by='Importance',ascending=False).head(25).loc[:,'Feature'].tolist()
                X_cols_classi_metrics=[col for col in X_cols_classi if col in X_cols_classi_metrics_selected]            
                grid_search_classi,data_preds=run_RF_classi(data_classi_train,X_cols_classi_metrics,y_coln_classi,
                        test_size=0.34,data_test=data_classi_test,data_out_fh=grid_search_classi_metrics_fh) 
                get_RF_classi_metrics(grid_search_classi_metrics_fh,data_dh='data_ml/',plot_dh='plots/')
                
            y_coln_classi="classes"
            feature_importances_classi_fh="%s_%s_.csv" % (grid_search_classi_fh,'featimps')
            feature_importances_classi=pd.read_csv(feature_importances_classi_fh)
            data_classi_test=pd.read_csv(data_fh.replace("data_ml_","data_ml_classi_test_"))
            data_classi_test  =data_classi_test.set_index("mutids",drop=True)
            data_classi_train=pd.read_csv(data_fh.replace("data_ml_","data_ml_classi_train_"))
            data_classi_train  =data_classi_train.set_index("mutids",drop=True)
            #regress
            y_coln_regress="FCA_norm"
            data_regress_train=data_classi_train                
            data_regress_train=X_cols2binary(data_regress_train,[y_coln_classi])
            data_fit=set_index(data_fit,"mutids")
            data_regress_train.loc[:,y_coln_regress]=\
            data_fit.loc[data_classi_train.index.values,y_coln_regress]
            data_regress_train=denan(data_regress_train,axis='both',condi='all any')
            data_regress_test=data_classi_test                
            if y_coln_classi in data_regress_test.columns.tolist(): 
                data_regress_test=data_regress_test.drop(y_coln_classi,axis=1)
            else:
                data_regress_test=denan(data_regress_test,axis='both',condi='all any')
                data_regress_test=X_cols2binary(data_regress_test,[y_coln_classi])
            data_regress_test=denan(data_regress_test,axis='both',condi='all any')
            X_cols_regress_class_fit=[]
            X_cols_regress_selected=feature_importances_classi\
            .sort_values(by='Importance',ascending=False).head(25).loc[:,'Feature'].tolist()
            X_cols_regress=[col for col in X_cols_regress_class_fit+X_cols_regress_selected\
                            if col in data_regress_test.columns.tolist()]
            if y_coln_regress in X_cols_regress:
                X_cols_regress.remove(y_coln_regress)
            data_regress_train.to_csv(data_fh.replace("data_ml_","data_ml_regress_train_"))
            data_regress_test.to_csv(data_fh.replace("data_ml_","data_ml_regress_test_"))
            # try:
            grid_search_regress_metrics,data_preds_regress_metrics=\
            run_RF_regress(data_regress_train,X_cols_regress,y_coln_regress,
                            test_size=0.34,data_test=data_regress_test,data_out_fh=grid_search_regress_metrics_fh)
            get_RF_regress_metrics(grid_search_regress_metrics_fh,data_dh='data_ml/',plot_dh='plots/')
            grid_search_regress,data_preds_regress=\
            run_RF_regress(data_regress_train,X_cols_regress,y_coln_regress,
                            test_size=0,data_test=data_regress_test,data_out_fh=grid_search_regress_fh)
            data_preds_regress.to_csv(data_fh.replace("data_ml_","data_ml_regress_preds_"))
            data_regress_all=data_preds_regress.append(data_regress_train)
            data_regress_all.to_csv(data_fh.replace("data_ml_","data_ml_regress_all_"))
            data_regress2data_fit(prj_dh,data_fit_key,data_regress_all)
            return data_regress_all
            # except:
                # logging.info("skipping: %s : requires more data" % basename(data_fit_key))
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
    
    
    