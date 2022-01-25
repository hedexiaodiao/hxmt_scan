#!/home/hxmt/nangyi/anaconda2/bin/python
from xspec import *
import math
import numpy as np
import time
from astropy.io import fits as pf
import sys
import matplotlib.pyplot as plt
plt.switch_backend('agg')

time1= time.time()

###########Class and Function defined by HDPC################################
sys.path.append("/hxmt/work/HXMT_scan_data/psfcrt/tools")#lib Wpath
from readxml import *

import createpha_beta as createpha
import psfME_gen as hxmtpsf

print "import class success"
if len(sys.argv)<2:
    print "Need the config file!"
    sys.exit(1)

print sys.argv[0],sys.argv[1]
cfg = sys.argv[1]
readcfg = loadDom(cfg)
infilestr = readcfg.getTagText("infilenamelist")
inpathstr = readcfg.getTagText("inpath")
outpathstr = readcfg.getTagText("outpath")
outfilestr = readcfg.getTagText("outfilename")
print inpathstr,infilestr,outpathstr
inpathstr = inpathstr.strip()
infilestr0 = infilestr.split()[0]
infilestr1 = infilestr.split()[1]
infilestr2 = infilestr.split()[2]
infilestr3 = infilestr.split()[3]
infile0 = (inpathstr+infilestr0).encode()
infile1 = (inpathstr+infilestr1).encode()
infile2 = (inpathstr+infilestr2).encode()
infile3 = (inpathstr+infilestr3).encode()
filename = filter(str.isdigit,infile3.split("/")[-1])
outpathstr = outpathstr.strip()
outfilestr0 = outfilestr.split()[0]
outfile0 = (outpathstr+outfilestr0+".log").encode()
outfile1 = (outpathstr+outfilestr0+".par").encode()
outfile2 = (outpathstr+outfilestr0+".eps").encode()
outfile3 = (outpathstr+outfilestr0).encode()
outfile4 = (outpathstr+outfilestr0+".src").encode()
outfile5 = (outpathstr+outfilestr0).encode()#.split("/")[-1]
print ("the count rate file : ",infile0,infile1,infile2)
print ("the out put file : ",outfile0,outfile1,outfile2)
##############load the hxmtpsf model#################################
try:
    logFile = Xset.openLog(outfile0)
    parFile = open(outfile1,"w")
except:
    print "Open log or par file failed!"
    sys.exit(1)
print infile0,infile1,infile2,hxmtpsf.infile
AllModels.addPyMod(hxmtpsf.hxmtpsf, hxmtpsf.myModelParInfo, 'add')
print "add hxmtpsf to the model"
#########calculate the src map in the region of interesting ################
AllModels.clear()

temp= "powerlaw+hxmtpsf"
k = 1

parnumb = (k-1)*5+7

#######################Generate fake pha for xspec##########################
createpha.cpha(infile0,"%s0"%outfile5)
createpha.cpha(infile1,"%s1"%outfile5)
createpha.cpha(infile2,"%s2"%outfile5)

###############load data and fit###########################################
AllData("1:1 %s0.grp 2:2 %s1.grp 3:3 %s2.grp"%(outfile5,outfile5,outfile5))
print "alldata loaded" ,infilestr
print "Fit Begin"

d1 = AllData(1)
d2 = AllData(2)
d3 = AllData(3)

AllModels+=temp

m1 = AllModels(1)
m2 = AllModels(2)
m3 = AllModels(3)
i=0
Fit.query = "yes"
m1.setPars({i*5+3:" 1,-1"})
m1.setPars({i*5+4:"60,-0.001"})
m1.setPars({i*5+5:"83.633,-1.0"})
m1.setPars({i*5+6:"22.015,-1.0"})
m1.setPars({i*5+7:"200,1.0"})

m2.setPars({1:"0 -1",2:"10 0.1"})
m2.setPars({i*5+4:"0 -0.001"})

m3.setPars({1:"0 -1",2:"10 0.1"})
m3.setPars({i*5+4:"-60  -0.001"})

print "****************** First max 50 FIT*********************"
print "****************** First max 50 FIT*********************"

Fit.perform()
Fit.error("maximum 150 1.0 1-%s"%(parnumb+1))
fig,axe = plt.subplots(3,sharex=True, sharey=False,figsize=(24,15))
axe[0].plot(d1.noticed, (np.array(d1.values)),'b-',lw=1,label='data')
axe[0].plot(d1.noticed, np.array(m1.folded(1)),'r-',lw=2,label='PSF')
lgd = axe[0].legend(loc='upper right',fontsize=10)
axe[0].set_title(filename,fontsize=12)
axe[1].plot(d2.noticed, (np.array(d2.values)),'b-',lw=1)
axe[1].plot(d2.noticed, (np.array(d2.values))-np.array(m2.folded(2))-50,'k-',lw=2)
axe[1].plot(d2.noticed, m2.folded(2),'r-',lw=2)
axe[1].set_ylabel('counts/sec',fontsize=20)
axe[2].plot(d3.noticed, (np.array(d3.values)),'b-',lw=1)
axe[2].plot(d3.noticed, (np.array(d3.values))-np.array(m3.folded(3))-50,'k-',lw=2)
axe[2].plot(d3.noticed, np.array(m3.folded(3)),'r-',lw=2)

plt.savefig('%s_LcFit.png'%(filename),dpi = 500)

m1(2).frozen = True
m2(2).frozen = True
m3(2).frozen = True

bkg1 = m1(2).values[0]
bkg2 = m2(2).values[0]
bkg3 = m3(2).values[0]
bkg = [bkg1,bkg2,bkg3]
fdtemp1 = np.array(m1.folded(1))
fdtemp2 = np.array(m2.folded(2))
fdtemp3 = np.array(m3.folded(3))

m1.setPars({7:"0 -1"})
#Fit.perform()

contri1 = fdtemp1 - np.array(m1.folded(1))
contri2 = fdtemp2 - np.array(m2.folded(2))
contri3 = fdtemp3 - np.array(m3.folded(3))
times1 = np.array(d1.noticed)[((fdtemp1 - np.array(m1.folded(1)))>0)]
times2 = np.array(d2.noticed)[((fdtemp2 - np.array(m2.folded(2)))>0)]
times3 = np.array(d3.noticed)[((fdtemp3 - np.array(m3.folded(3)))>0)]

srcpot=[]

try:
    a=((times1[1:]-times1[:-1])!=1)
    c=np.nonzero(a)[0]
    pot1=np.r_[0,c+1]
    pot2=np.r_[c,times1.shape[0]-1]
    pot1 = times1[pot1]
    pot2 = times1[pot2]
    pot=[pot1,pot2]
    srcpot.append(pot)
except:
    srcpot.append([])

try:
    a=((times2[1:]-times2[:-1])!=1)
    c=np.nonzero(a)[0]
    pot1=np.r_[0,c+1]
    pot2=np.r_[c,times2.shape[0]-1]
    pot1 = times2[pot1]
    pot2 = times2[pot2]
    pot=[pot1,pot2]
    srcpot.append(pot)
except:
    srcpot.append([])

try:
    a=((times3[1:]-times3[:-1])!=1)
    c=np.nonzero(a)[0]
    pot1=np.r_[0,c+1]
    pot2=np.r_[c,times3.shape[0]-1]
    pot1 = times3[pot1]
    pot2 = times3[pot2]
    pot=[pot1,pot2]
    srcpot.append(pot)
except:
    srcpot.append([])

AllData.clear()

Fit.query = "yes"
#AllModels+=temp
listt = [0,1,2]
rollt = [60,0,-60]
fragnum=[6,26,6]
for lcnum in listt:
    AllData("%s%s.grp"%(outfile5,lcnum))
    NormLc=[]
    AlfaLc=[]
    erroLc=[]
    if len(srcpot[lcnum])==0:
        np.save('Angle_normT_%s.npy'%(lcnum),np.array([[],[],[]]))
        continue
    for pot in xrange(len(srcpot[lcnum][0])):
        start = srcpot[lcnum][0][pot]
        stop = srcpot[lcnum][1][pot]
        if stop-start < 1:
            continue
        AllData.ignore("*")
        AllData.notice("%i-%i"%(start,stop))
        m1 = AllModels(1)
        d1 = AllData(1)
        i=0
        m1.setPars({1:"0 -1",2:"%s -0.1"%bkg[lcnum]})
        Fit.query = "yes"
        m1.setPars({i*5+3:" 1   -1"})
        m1.setPars({i*5+4:"%s    -0.001"%(rollt[lcnum])})
        m1.setPars({i*5+5:"83.633,-1.0"})
        m1.setPars({i*5+6:"22.015,-1.0"})
        m1.setPars({i*5+7:"100,1.0,0.,0.,100000.,100000,"})
        Fit.perform()
        Fit.error("maximum 100 1.0 1-7")
        Fit.perform()
        alpha = np.array(hxmtpsf.D_ra)[start-1:stop+1]
        beta = np.array(hxmtpsf.D_dec)[start-1:stop+1]
        alpha0 = np.array(hxmtpsf.D_dec0)[start-1:stop+1]
        print len(beta),len(alpha0)
        step = np.linspace(-1.,1.,fragnum[lcnum]+1)
        if lcnum==1:
            step = np.linspace(-4.,4.,fragnum[lcnum]+1)
        for j in xrange(fragnum[lcnum]):
            seq = (alpha>step[j])&(alpha<step[j+1])&(beta>-4)&(beta<4)
            if lcnum==1:
                seq = (beta>step[j])&(beta<step[j+1])&(alpha>-1.)&(alpha<1.)
            a=np.nonzero(seq)[0]
            if len(a)>3:
                startT = start + a[0]
                stopT = start + a[-1]
                #print stop,start,a[0],a[-1]
                AllData.ignore("*")
                AllData.notice("%i-%i"%(startT,stopT))
                Fit.perform()
                Fit.error("maximum 100 1.0 1-7")
                alphaT = hxmtpsf.D_dec0[startT-1:stopT] if lcnum!=1 else hxmtpsf.D_ra0[startT-1:stopT]
                betaT = hxmtpsf.D_dec[startT-1:stopT]
                errT = (m1(7).error[1]-m1(7).error[0])/2.
                if errT==0:
                    print "ERROR::: err==0!!!"
                    print "Norm:",m1(7).values,m1(7).error
                    with open('Log.txt','w') as f:
                        print>>f,"Norm:",m1(7).values,'\t',m1(7).error
                    float('l')
                erroLc.append(errT)
                NormLc.append(m1(7).values[0])
                AlfaLc.append(np.mean(alphaT))
            else:
                erroLc.append(-99)
                NormLc.append(-99)
                AlfaLc.append(-99.)
            #print alphaT,betaT
            #print m1(7).values[0]
            #float(raw_input("::"))
    np.save('Angle_normT_%s.npy'%(lcnum),(AlfaLc,NormLc,erroLc))
    AllData.clear()

parFile.close()
time2 = time.time()
print 'fit finished ',time2
print "time expended: ",(time2-time1)/3600.0,' hours'

