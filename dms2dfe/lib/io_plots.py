#!usr/bin/python

# Copyright 2016, Rohan Dandage <rraadd_8@hotmail.com,rohan@igib.in>
# This program is distributed under General Public License v. 3.  

import matplotlib.pyplot as plt
import numpy as np

def saveplot(plot_fh,form='both',
            transparent=True,
            tight_layout=True,):
    if not plot_fh is None:
        def save(plot_fh,form,transparent):
            if '.%s' % form in plot_fh:
                plot_out_fh=plot_fh
            else:
                plot_out_fh='%s.%s' % (plot_fh,form)
            print plot_out_fh
            plt.savefig(plot_out_fh,format=form,transparent=transparent)
        if tight_layout:
            plt.tight_layout()
        if plot_fh!=None:
            if form=="pdf" or form=="both":
                save(plot_fh,'pdf',transparent)
            if form=="png" or form=="both":
                save(plot_fh,'png',transparent)

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