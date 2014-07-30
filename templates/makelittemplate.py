#!/usr/bin/env python

import sys,os,glob
#,re,numpy,math,pyfits,glob,shutil,glob
import optparse
import numpy as np
import pylab as pl
from scipy import stats as stats 
from scipy import optimize
import math 
from mpmath import polyroots
import time
from scipy.signal import cspline1d, cspline1d_eval
import pickle as pkl


from templutils import *

allbands = ['U','B','V','R','I','g','r','i','z']



if __name__=='__main__':
    parser = optparse.OptionParser(usage="makelittemplate.py ", conflict_handler="resolve")
    parser.add_option('--showonly', default=False, action="store_true",
                      help='show only')
    parser.add_option('-b','--band', default='V' , type="string",
                      help='degree of polinomial fit')
    parser.add_option('-i','--inst', default='all' , type="string",
                      help='instrument to be used: shooter, mini, keploer, or all')
    options,  args = parser.parse_args()
    if len(args)>1:
        sys.argv.append('--help')
    
        options,  args = parser.parse_args()
        sys.exit(0)
        
    if options.band == 'V': band = 'V'
    else:band = options.band
    if band not in allbands:
        print "band "+band+" not available"
        sys.exit()

        
    sne,templates=loadlitlist(band)
    sne=splinetemplates(sne)
    ax=pl.figure()
    legends=[]

    
    for s in sne:         
        pl.errorbar(s.phot[0]-s.mjdmax,s.normphot,s.phot[2],fmt='ko')
        leg,=pl.plot(s.new_x,s.new_y)
#        pl.show()
        legends.append(leg)
        
   # print ([s['new_y'] for s in sne])
    mytemplate=open("mytemplate"+band+".dat",'w')        
    mt=Mytempclass()
    pklfile = open('mytemplate'+band+'.pkl', 'wb')
    mt.loadtemplate(band, x=sne[0].new_x,mean=smoothListGaussian(stats.stats.nanmean([s.new_y for s in sne],axis=0)),median=smoothListGaussian(stats.stats.nanmedian([s.new_y for s in sne],axis=0)),std=smoothListGaussian(stats.stats.nanstd([s.new_y for s in sne],axis=0)))
    mt.template[band].mean=mt.template[band].mean[18:]
    mt.template[band].median=mt.template[band].median[18:]
    mt.template[band].std=mt.template[band].std[18:]
    mt.template[band].x=mt.template[band].x[18:]
    pkl.dump(mt.template[band],pklfile)

#    pl.fill_between(mt.template[band].x[100:],mt.template[band].mean[100:]-mt.template[band].std[100:],mt.template[band].mean[100:]+mt.template[band].std[100:],facecolor='gray', alpha=0.5)#, transform=trans)
    pl.plot ( mt.template[band].x,mt.template[band].mean+3,'k-',linewidth=3)
    pl.plot ( mt.template[band].x,mt.template[band].median+3,'b-',linewidth=3)
    pl.plot ( mt.template[band].x,mt.template[band].mean-mt.template[band].std+3,'b-',linewidth=1)
    pl.plot ( mt.template[band].x,mt.template[band].mean+mt.template[band].std+3,'b-',linewidth=1)

    pl.legend(legends,templates['sn'], loc=1, ncol=1,prop={'size':8})        
    pl.ylim(7,-1)
    pl.xlim(-12,42)
#    pl.show()

    pl.savefig("template"+band+".png")

    print >> mytemplate, "#phase, mean, median, stdev"
    for i,x in enumerate(mt.template[band].x):
        print >> mytemplate, x,mt.template[band].mean[i],mt.template[band].median[i],mt.template[band].std[i]

    
    
        
