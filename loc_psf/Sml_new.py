#!/home/hxmt/nangyi/anaconda2/bin/python
from xspec import *
import math
import numpy as np
import time
from astropy.io import fits as pf
import sys,os
import matplotlib.pyplot as plt
#plt.switch_backend('agg')

time1= time.time()

###########Class and Function defined by HDPC################################
#sys.path.append("/sharefs/hbkg/data/SCAN/PSF/psf_201904")#lib path
#from readxml import *
#import createpha_beta as createpha
#import hxmt_paraodic as hxmtpsf
sys.path.append("/sharefs/hbkg/user/nangyi/lick_tools")#lib path
from readxml import *
#import createpha_beta as createpha
import hxmt_paraodic3d as hxmtpsf

print "import class success"
#######read input output template paths and files from  config file##########
def findpot(inpu):
    times = np.copy(inpu)
    try:
        a=((times[1:]-times[:-1])!=1)
        c=np.nonzero(a)[0]
        pot1=np.r_[0,c+1]
        pot2=np.r_[c,times.shape[0]-1]
        pot1 = times[pot1]
        pot2 = times[pot2]
        pot=[pot1,pot2]
        return pot
    except:
        return []

########################################################################
if len(sys.argv)<2:
    print "Need the config file!"
    sys.exit(1)

start = int(sys.argv[-1])
step = int(sys.argv[-2])
paths = '/sharefs/hbkg/data/SCAN/PSF/HE/crabtxt/txt_para/'
paths = '/sharefs/hbkg/data/SCAN/PSF/HE/crabtxt2020/2018/'
paths = '/hxmt/work/HXMT_scan_data/psfcrt/LE/crabtxt/2018/'
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
filename = '0101'+os.getcwd().split("/")[-1] #filter(str.isdigit,infile3.split("/")[-1])
#filename = filter(str.isdigit,infile3.split("/")[-1])
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

temp= "constant*hxmtpsf"
k = 1

parnumb = (k-1)*5+7
#######################Generate fake pha for xspec##########################
#createpha.cpha("hepoint0.fits","hepoint0")
#createpha.cpha("hepoint1.fits","hepoint1")
#createpha.cpha("hepoint2.fits","hepoint2")
#createpha.cpha("mepoint0.fits","mepoint0")
#createpha.cpha("mepoint1.fits","mepoint1")
#createpha.cpha("mepoint2.fits","mepoint2")
#createpha.cpha("lepoint0.fits","lepoint0")
#createpha.cpha("lepoint1.fits","lepoint1")
#createpha.cpha("lepoint2.fits","lepoint2")


###############load data and fit###########################################
#AllData("1:1 %s0.grp 2:2 %s1.grp 3:3 %s2.grp"%(outfile5,outfile5,outfile5))
AllData("1:1 hepoint0.grp 2:2 hepoint1.grp 3:3 hepoint2.grp 4:4 mepoint0.grp 5:5 mepoint1.grp 6:6 mepoint2.grp 7:7 lepoint0.grp 8:8 lepoint1.grp 9:9 lepoint2.grp")

lth = len(AllData(1).values)
print "alldata loaded" ,infilestr
print "Fit Begin"

d1 = AllData(1)
d2 = AllData(2)
d3 = AllData(3)
d4 = AllData(4)
d5 = AllData(5)
d6 = AllData(6)
d7 = AllData(7)
d8 = AllData(8)
d9 = AllData(9)

AllModels+=temp

m1 = AllModels(1)
m2 = AllModels(2)
m3 = AllModels(3)
m4 = AllModels(4)
m5 = AllModels(5)
m6 = AllModels(6)
m7 = AllModels(7)
m8 = AllModels(8)
m9 = AllModels(9)

Fit.query = "yes"
i=0
ra,dec = 83.6330,22.015
m1.setPars({1:"1 -1",2:"0 -0.1"})
m2.setPars({1:"1 -1",2:"0 -0.1"})
m3.setPars({1:"1 -1",2:"0 -0.1"})

m1.setPars({i*5-1+4:"60,-0.001"})
m1.setPars({i*5-1+5:"%s,-1.0"%ra})
m1.setPars({i*5-1+6:"%s,-1.0"%dec})
m2.setPars({i*5-1+4:"0,-0.001"})
m3.setPars({i*5-1+4:"-60,-0.001"})
m2.setPars({i*5-1+3:" 0,-1"})
m3.setPars({i*5-1+3:" 0,-1"})
m1.setPars({6:"150., 0.1"})

m4.setPars({1:"1 -1",2:"1, -0.1"})
m5.setPars({1:"1 -1",2:"1, -0.1"})
m6.setPars({1:"1 -1",2:"1, -0.1"})

m4.setPars({i*5-1+4:"60,-0.001"})
m5.setPars({i*5-1+4:"0,-0.001"})
m6.setPars({i*5-1+4:"-60,-0.001"})
m4.setPars({6:"127., 0.1"})
m5(6).link='24'
m6(6).link='24'

m7.setPars({1:"1 -1",2:"2  -0.1"})
m8.setPars({1:"1 -1",2:"2  -0.1"})
m9.setPars({1:"1 -1",2:"2  -0.1"})

m7.setPars({i*5-1+4:"60,-0.001"})
m8.setPars({i*5-1+4:"0,-0.001"})
m9.setPars({i*5-1+4:"-60,-0.001"})
m7.setPars({6:"130., 0.1"})
m8(6).link='42'
m9(6).link='42'

m1(4).values=[83.633,0.01,ra-1.,ra-1.,ra+1,ra+1]
m1(5).values=[22.014,0.01,dec-1,dec-1,dec+1,dec+1]

print "****************** First max 50 FIT*********************"
print "****************** First max 50 FIT*********************"

AllData.ignore('**-29000')
AllData.notice('%s-%s'%(start,start+step-1))
AllData.ignore('29000-**')
#Fit.perform()
Fit.steppar('4 83.533 83.733 20 5 21.915 22.115 20')
stat=Fit.stepparResults('statistic');gra = Fit.stepparResults('4');gdec = Fit.stepparResults('5')
idx = np.argmin(stat)
m1(4).values=[gra[idx],-0.01];m1(5).values=[gdec[idx],-0.01]
Fit.perform()
m1(4).values=[gra[idx],0.01];m1(5).values=[gdec[idx],0.01]
Fit.perform()
if (m1(4).values[0]==83.533)|(m1(4).values[0]==83.733)|(m1(5).values[0]==21.915)|(m1(5).values[0]==22.115):
    Fit.steppar('4 83.233 84.033 20 5 21.615 22.415 20')
    Fit.perform()
Fit.error("1.0 1-54")
ra = m1(4).values[0]
dec = m1(5).values[0]
with open('pointinfonew.txt','a')as f:
    print>>f,step,start,ra,m1(4).error[0],m1(4).error[1],dec,m1(5).error[0],m1(5).error[1],Fit.statistic,Fit.dof
#Fit.perform()

with open('poinnew.log','a')as f:
    print>>f,step,start

