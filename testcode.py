from snclasses import *
from templutils import *
import optparse
import readinfofile as ri
import pylabsetup
su=setupvars()

f=glob.glob(os.environ['SESNPATH']+"/finalphot/*05bf*[cf]")
assert len(f)>0, "no such SN in the data archive"
f=f[0]
##note: all files have extention .c or .f in cfa3 and cfa4

Vmax,sntype=ri.readinfofile()
thissn=mysn(f)
fnir = True		
lc,flux,dflux, snname = thissn.loadsn(f,fnir, verbose=True)
assert  snname.lower() in Vmax.keys(),  "skipping object missing from metadata table %s"%snname.lower()
     
thissn.setsn(sntype[snname.lower()],Vmax[snname.lower()])
thissn.setphot()
##look for galaxy extinction correction
myebmv=0
for snoff in su.ebmvs.iterkeys():
    if thissn.name.endswith((snoff.strip()).lower()):
            myebmv=su.ebmvs[snoff]

thissn.getphot(myebmv)


thissn.getcolors()
thissn.printsn(photometry=True, fout="sn05bf.phot.lcv")
thissn.printsn(color=True, fout="sn05bf_color.phot.lcv")

thissn.plotsn(photometry=True,  offsets=True, aspect=0.5, Vmax=True, savepng=True, show=True,  legendloc=2)
thissn.plotsn(color=True,show=True,  offsets=True, aspect=0.5, Vmax=True, savepng=True)
