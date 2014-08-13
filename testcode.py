from snclasses import *
from templutils import *
import optparse
import readinfofile as ri
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
for snoff in ebmvs.iterkeys():
    if thissn.name.endswith((snoff.strip()).lower()):
            myebmv=ebmvs[snoff]


thissn.printsn_textable(photometry=True, fout=snname+".phot.tex")   
thissn.getphot(myebmv)
thissn.getcolors()
thissn.printsn(photometry=True)
thissn.printsn(color=True)
thissn.plotsn(photometry=True,show=False, fig=0,  relim=False, offsets=True, aspect=0.5, Vmax=False, save=True)
thissn.plotsn(color=True,show=False, fig=1, ylim=(maxmag,minmag), xlim=(mint-10,maxt+10), relim=False, offsets=True, mylabel=ylabel, aspect=0.5, Vmax=False, legendloc=legendloc,noylabel=noylabel, save=True)
