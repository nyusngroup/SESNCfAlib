#!/usr/bin/env python
import sys,os,glob, inspect
#,re,numpy,math,pyfits,glob,shutil,glob
import optparse
import scipy as sp
import numpy as np
from scipy import optimize
from scipy.interpolate import interp1d
from mpmath import polyroots
import time
cmd_folder = os.path.realpath(os.path.abspath(os.path.split(inspect.getfile( inspect.currentframe() ))[0]) + "/templates")
if cmd_folder not in sys.path:
     sys.path.insert(0, cmd_folder)
import copy
import pickle as pkl
from snclasses import *
from templutils import *
from utils import *
from fitutils import *
from plotutils import *
from random import shuffle
from fitgauss2sntemplate import mygauss,exprise
cmd_folder = os.path.realpath(os.path.abspath(os.path.split(inspect.getfile( inspect.currentframe() ))[0]) + "/getCfASNtable")
if cmd_folder not in sys.path:
     sys.path.insert(0, cmd_folder)
from getfromtable import grepstuff,makeintodictionary
#pl.ion()
#LN10x2p5=5.75646273249

from cosmdist import cosmo_dist

info= grepstuff ('fedhere','nel1962','z',None)#sn.strip())
#  print  info
zdic= makeintodictionary(info)
print zdic
sne= zdic.keys()

for i,s in enumerate(sne):
    print s,zdic[s]['z'],
    try:
        print cosmo_dist([0],[float(zdic[s]['z'])],lum=1,Mpc=1)[0]
    except:
        pass
