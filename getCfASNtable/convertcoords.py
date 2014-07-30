from getfromtable import grepstuff
from myastrotools import raseg2deg,decseg2deg

#list = open("../snlist").readlines()

#for sn in list:
info= grepstuff ('fedhere','nel1962','ra,dec',None)#sn.strip())
print info[1][0]
for i in range(len(info[1][0])):
    print info[1][0][i],
    print   raseg2deg(info[1][1][i]) ,
    print decseg2deg(info[1][2][i])

    
