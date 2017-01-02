#!usr/bin/python

# Copyright 2016, Rohan Dandage <rraadd_8@hotmail.com,rohan@igib.in>
# This program is distributed under General Public License v. 3.  

import matplotlib
matplotlib.style.use('ggplot')
matplotlib.rcParams['axes.unicode_minus']=False
matplotlib.use('Agg') # no Xwindows
import matplotlib.pyplot as plt

import seaborn as sns
import numpy as np

from dms2dfe.lib.io_dfs import debad

def saveplot(plot_fh,form='both',
            transparent=False,
            tight_layout=True,
            print_fh=False):
    if not plot_fh is None:
        def save(plot_fh,form,transparent):
            if '.%s' % form in plot_fh:
                plot_out_fh=plot_fh
            else:
                plot_out_fh='%s.%s' % (plot_fh,form)
            if print_fh:
                print plot_out_fh
            plt.savefig(plot_out_fh,format=form,transparent=transparent)
        if tight_layout:
            plt.tight_layout()
        if plot_fh!=None:
            if form=="pdf" or form=="both":
                save(plot_fh,'pdf',transparent)
            if form=="png" or form=="both":
                save(plot_fh,'png',transparent)
        plt.clf()
        plt.close()

def get_axlims(X,Y,space=0.2,equal=False):
    xmin=np.min(X)
    xmax=np.max(X)
    xlen=xmax-xmin
    ymin=np.min(Y)
    ymax=np.max(Y)
    ylen=ymax-ymin
    xlim=(xmin-space*xlen,xmax+space*xlen)
    ylim=(ymin-space*ylen,ymax+space*ylen)
    if not equal:
        return xlim,ylim
    else:
        lim=[np.min([xlim[0],ylim[0]]),np.max([xlim[1],ylim[1]])]
        return lim,lim

from scipy import stats
def plot_scatter_reg(data_all,cols,
                     xlabel=None,ylabel=None,title=None,
                     color_scatter="gray",color_line='k',
                     logscale=False,
                     space=0.2,figsize=[2,2],
                     ax=None,plot_fh=None):
    if ax==None:
        fig=plt.figure(figsize=figsize,dpi=300)
        ax=plt.subplot(111)
    data_all=debad(data_all.loc[:,cols],axis=0,condi='any',bad='nan')
    if logscale:
        data_all=debad(data_all,axis=0,condi='any',bad=0)
        data_all=data_all.apply(np.log2)    
    ax=sns.regplot(data=data_all,x=cols[0],y=cols[1],
                line_kws={"color":color_line},
                scatter_kws={"color":color_scatter},
                ax=ax)
    r, _ = stats.pearsonr(data_all.loc[:,cols[0]],
                         data_all.loc[:,cols[1]])
    ax.legend(["r=%.2f" % (r)],loc="upper left")
    if not title is None:
        sns.plt.title(title)
    xlim,ylim=get_axlims(data_all.loc[:,cols[0]],
                 data_all.loc[:,cols[1]],space=space)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    saveplot(plot_fh,form='both',transparent=False)
    return ax
