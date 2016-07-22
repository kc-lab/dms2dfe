#!usr/bin/python

# Copyright 2016, Rohan Dandage <rraadd_8@hotmail.com,rohan@igb.in>
# This program is distributed under General Public License v. 3.  

import sys
from os.path import exists,splitext,basename
from os import makedirs
from glob import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from multiprocessing import Pool
import warnings
warnings.simplefilter(action = "ignore", category = FutureWarning)
warnings.simplefilter(action = "ignore", category = DeprecationWarning)
import logging
logging.basicConfig(format='[%(asctime)s] %(levelname)s\tfrom %(filename)s in %(funcName)s(..):%(lineno)d: %(message)s',level=logging.DEBUG) # filename=cfg_xls_fh+'.log'
from dms2dfe import configure
from dms2dfe.lib.io_ml import data_fit_feats2combo,y2classes,X_cols2numeric,denanrows,run_RF,data_fit2ml
from dms2dfe.lib.global_vars import mut_types_form

def main(prj_dh):
    """
    This modules trains a Random Forest classifier with given data in `data_fit` format.
    
    This plots the results in following visualisations.
    
    .. code-block:: text
    
        ROC plots
        Relative importances of features
    
    :param prj_dh: path to project directory.
    """    
    logging.info("start")

    if not exists(prj_dh) :
        logging.error("Could not find '%s'" % prj_dh)
        sys.exit()
    configure.main(prj_dh)
    from dms2dfe.tmp import info
    cores=int(info.cores)
    
    global prj_dh_global,data_feats,y_coln
    prj_dh_global=prj_dh
    type_form="aas"
    if not exists("%s/plots/%s" % (prj_dh,type_form)):
        makedirs("%s/plots/%s" % (prj_dh,type_form))
    if not exists("%s/data_ml/%s" % (prj_dh,type_form)):
        makedirs("%s/data_ml/%s" % (prj_dh,type_form))

    data_feats=pd.read_csv("%s/data_feats/aas/feats_all" % (prj_dh))
    if "aasi" in data_feats.columns:
        data_feats=data_feats.set_index("aasi",drop=True)
    y_coln='class_fit'
    data_fit_keys = ["data_fit/%s/%s" % (type_form,basename(fh)) \
                     for fh in glob("%s/data_fit/aas/*" % prj_dh) \
                     if not "infered" in basename(fh)]
    data_fit_keys = np.unique(data_fit_keys)
    if len(data_fit_keys)!=0:
        # for data_fit_key in data_fit_keys:
        #     pooled_io_ml(data_fit_key)
        # pooled_io_ml(data_fit_keys[-1])
        pool_io_ml=Pool(processes=int(cores)) 
        pool_io_ml.map(pooled_io_ml,data_fit_keys)
        pool_io_ml.close(); pool_io_ml.join()
    else:
        logging.info("already processed")
#     y_coln='class_comparison'
#     data_fit_keys = ["data_comparison/%s/%s" % (type_form,basename(fh)) for fh in glob("%s/data_comparison/aas/*" % prj_dh)]
#     data_fit_keys = np.unique(data_fit_keys)
#     if len(data_fit_keys)!=0:
# #         pooled_io_ml(data_fit_keys[0])
#         pool_io_ml=Pool(processes=int(cores)) 
#         pool_io_ml.map(pooled_io_ml,data_fit_keys)
#         pool_io_ml.close(); pool_io_ml.join()
#     else:
#         logging.info("already processed")
    logging.shutdown()

def pooled_io_ml(data_fit_key):
    """
    This module makes use of muti threading to speed up `dms2dfe.lib.io_ml.data_fit2ml`.     
    
    :param data_fit_key: in the form <data_fit>/<aas/cds>/<name of file>.
    """
    data_fit2ml(data_fit_key,prj_dh_global,data_feats)
    
if __name__ == '__main__':
    main(sys.argv[1])