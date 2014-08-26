import numpy as np
import glob, pickle
import os,inspect,sys

try:
     libpath=os.environ['SESNPATH']

except KeyError:
     print "must set environmental variable SESNPATH"
     sys.exit()

cmd_folder = os.path.realpath(libpath+"/SESNCFAlib")
if cmd_folder not in sys.path:
     sys.path.insert(0, cmd_folder)


cmd_folder =  os.path.realpath(libpath+"/SESNCFAlib/templates")
if cmd_folder not in sys.path:
     sys.path.insert(0, cmd_folder)


from snclasses import *
from templutils import *
import optparse
import pylabsetup

su=setupvars()

f=glob.glob(os.environ['SESNPATH']+"/finalphot/*05bf*[cf]")
assert len(f)>0, "no such SN in the data archive"
f=f[0]
##note: all files have extention .c or .f in cfa3 and cfa4

thissn=mysn(f)
fnir = True		
lc,flux,dflux, snname = thissn.loadsn(f,fnir, verbose=True)
#assert  snname.lower() in Vmax.keys(),  "skipping object missing from metadata table %s"%snname.lower()
thissn.readinfofileall(verbose=False, earliest=False, loose=True)
     
thissn.setsn(thissn.metadata['Type'],thissn.Vmax)
thissn.setphot()
##look for galaxy extinction correction
myebmv=0

assert  thisebmv=su.ebmvs[thissn.snnameshort], "couldnt find ebmv correction in 
    su=setupvars()"
try:
     thisebmv+=su.ebmvhost[thissn.snnameshort]
except KeyError:
     pass 
#for snoff in su.ebmvs.iterkeys():
#    if thissn.name.endswith((snoff.strip()).lower()):
#            myebmv=su.ebmvs[snoff]

thissn.getphot(myebmv)


thissn.getcolors()
thissn.printsn(photometry=True, fout="sn05bf.phot.lcv")
thissn.printsn(color=True, fout="sn05bf_color.phot.lcv")

thissn.plotsn(photometry=True,  offsets=True, aspect=0.5, Vmax=True, savepng=True, show=True,  legendloc=2)
thissn.plotsn(color=True,show=True,  offsets=True, aspect=0.5, Vmax=True, savepng=True)
