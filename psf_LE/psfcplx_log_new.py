#!/home/hxmt/nangyi/anaconda2/bin/python
from xspec import *
#import math
import numpy as np
import time
#from astropy.io import fits as pf
import sys,os
import matplotlib.pyplot as plt
#plt.switch_backend('agg')

time1= time.time()

###########Class and Function defined by HDPC################################

###sys.path.append("/hxmt/work/HXMT_scan_data/psfcrt/tools")#lib path
sys.path.append('/sharefs/hbkg/user/luoqi/GRB/work/ihep4/HXMT_GPS_Program/HXMT_GPS_py3Tools')
###sys.path.append(os.path.abspath(os.path.dirname(__file__))+'/HXMT_GPS_py3Tools')
from readxml import *
import createpha_beta as createpha
import cald_psf_sl as hxmtpsf
print("import class success")
#######read input output template paths and files from  config file##########

if len(sys.argv)<2:
    print("Need the config file!")
    sys.exit(1)

print(sys.argv[0],sys.argv[1])
cfg = sys.argv[1]
readcfg = loadDom(cfg)
infilestr = readcfg.getTagText("infilenamelist")
inpathstr = readcfg.getTagText("inpath")
outpathstr = readcfg.getTagText("outpath")
outfilestr = readcfg.getTagText("outfilename")
print(inpathstr,infilestr,outpathstr)
inpathstr = inpathstr.strip()
infilestr0 = infilestr.split()[0]
infilestr1 = infilestr.split()[1]
infilestr2 = infilestr.split()[2]
infilestr3 = infilestr.split()[3]
infile0 = (inpathstr+infilestr0)###.encode()
infile1 = (inpathstr+infilestr1)###.encode()
infile2 = (inpathstr+infilestr2)###.encode()
infile3 = (inpathstr+infilestr3)###.encode()
filename = filter(str.isdigit,infile3.split("/")[-1])
outpathstr = outpathstr.strip()
outfilestr0 = outfilestr.split()[0]
outfile0 = (outpathstr+outfilestr0+".log")###.encode()
outfile1 = (outpathstr+outfilestr0+".par")###.encode()
outfile2 = (outpathstr+outfilestr0+".eps")###.encode()
outfile3 = (outpathstr+outfilestr0).encode()
outfile4 = (outpathstr+outfilestr0+".src")###.encode()
outfile5 = (outpathstr+outfilestr0)###.encode()#.split("/")[-1]
print ("the count rate file : ",infile0,infile1,infile2)
print ("the out put file : ",outfile0,outfile1,outfile2)

outnpy_name = outpathstr+'/Angle_normC'
outpara_name = outpathstr+'/Para_log.txt'
print("The result will be generate into:\n",outnpy_name,outpara_name)
##############load the hxmtpsf model#################################
try:
    logFile = Xset.openLog(outfile0)
    parFile = open(outfile1,"w")
except:
    print("Open log or par file failed!")
    sys.exit(1)
print(infile0,infile1,infile2,hxmtpsf.infile)
AllModels.addPyMod(hxmtpsf.hxmtpsf, hxmtpsf.myModelParInfo, 'add')
print("add hxmtpsf to the model")
#########calculate the src map in the region of interesting ################
AllModels.clear()

temp= "powerlaw+hxmtpsf"
k = 1

parnumb = (k-1)*5+7
#######################Generate fake pha for xspec##########################
createpha.cpha(infile0,"%s0"%outfile5)
createpha.cpha(infile1,"%s1"%outfile5)
createpha.cpha(infile2,"%s2"%outfile5)
'''
#float('l')
###############load data and fit###########################################
AllData("1:1 %s0.grp 2:2 %s1.grp 3:3 %s2.grp"%(outfile5,outfile5,outfile5))
print "alldata loaded" ,infilestr
print "Fit Begin"
#AllData.ignore('81494-81700')
#AllData.ignore('75956-76200')
#AllData.ignore('50720-51300')
#AllData.ignore('64080-64120')
d1 = AllData(1)
d2 = AllData(2)
d3 = AllData(3)

AllModels+=temp

m1 = AllModels(1)
m2 = AllModels(2)
m3 = AllModels(3)

Fit.query = "yes"
i=0
m1.setPars({1:"0 -1",2:"10 0.1"})
m2.setPars({1:"0 -1",2:"10 0.1"})
m3.setPars({1:"0 -1",2:"10 0.1"})
m1.setPars({i*5+3:" 1,-1"})
m1.setPars({i*5+4:"60,-0.001"})
m1.setPars({i*5+5:"83.633,-1.0"})
m1.setPars({i*5+6:"22.015,-1.0"})
#m1.setPars({i*5+7:"0,-1.0"})
#m1.setPars({i*5+8:"0,-1.0"})
#m1.setPars({i*5+9:"0,-1.0"})

m2.setPars({i*5+4:"0 -0.001"})

m3.setPars({i*5+4:"-60  -0.001"})

print "****************** First max 50 FIT*********************"
print "****************** First max 50 FIT*********************"

Fit.perform()
#Fit.error("maximum 150 1.0 1-%s"%(parnumb+1))
fig,axe = plt.subplots(3,sharex=True, sharey=False,figsize=(24,15))
axe[0].plot(d1.noticed, (np.array(d1.values)),'b-',lw=1,label='data')
axe[0].plot(d1.noticed, (np.array(d1.values))-np.array(m1.folded(1))-50,'k-',lw=2)
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
plt.show()
float('l')

plt.savefig('%s_LcFit.png'%(filename),dpi = 500)

'''
#pqpqpqpqpqpqpqpqpqpqpqpqpqppqpqpqpqpqpqpqpqpqpqpqpqpqpqppqpqpqpqpqpqpqpqpqpqpqpqpqpqppqpqpqpqpqpqpqp
'''
fig,axe = plt.subplots(2,sharex=True, sharey=False,figsize=(15,10))
axe[0].plot(d1.noticed, (np.array(d1.values)),'b-',lw=1,label='data')
axe[0].plot(d1.noticed, np.array(m1.folded(1)),'r-',lw=2,label='PSF')
axe[0].grid()
axe[1].plot(d1.noticed, (np.array(d1.values))-np.array(m1.folded(1)),'k-',lw=2,label='Resdual')
axe[1].grid()
plt.show()

'''

AllModels.clear()
AllData.clear()

#pqpqpqpqpqpqpqpqpqpqpqpqpqppqpqpqpqpqpqpqpqpqpqpqpqpqpqppqpqpqpqpqpqpqpqpqpqpqpqpqpqppqpqpqpqpqpqpqp
Fit.query = "yes"
#AllModels+=temp
listt = [0,1,2]
rollt = [60,0,-60]
fragnum = [8,6,8]
psai = []
thta = []
phi = []
nmabs = []
i=0
#float('l')

os.system('rm %s*;rm %s;'%(outnpy_name,outpara_name))
crab_flux = 282

for lcnum in listt:
    AllData("%s%s.grp"%(outfile5,lcnum))
    AllModels+=temp
    m1 = AllModels(1)
    d1 = AllData(1)
    m1.setPars({1:"0 -1",2:"%s 0.1"%10})
    Fit.query = "yes"
    m1.setPars({i*5+3:" 1   -1"})
    m1.setPars({i*5+4:"%s    -0.001"%(rollt[lcnum])})
    m1.setPars({i*5+5:"83.633,-1.0"})
    m1.setPars({i*5+6:"22.015,-1.0"})
    m1.setPars({i*5+7:"0.,-0.01"})
    m1.setPars({i*5+8:"0.,-0.01"})
    m1.setPars({i*5+9:"0.,-0.01"})
    #m1.setPars({i*5+10:"200,1.0,0.,0.,1000.,1000,"})
    Fit.perform()
    m1.setPars({i*5+7:"0.,0.1"})
    m1.setPars({i*5+8:"0.,0.1"})
    m1.setPars({i*5+9:"0.,0.1"})
    Fit.perform()
    Fit.error("maximum 150 3.0 1-16")
    with open(outpara_name,'a')as f:
        for pa in range(1,17):
            try:
                print(pa,m1(pa).name,m1(pa).values[0],m1(pa).error[0],m1(pa).error[1],file=f)
            except:
                print(pa,m1(pa).name,m1(pa).values[0],-1,-1,file=f)
        print('*','CHIsqure',Fit.statistic,Fit.dof,file=f)
    m1.setPars({i*5+12:"0.,0.1"})
    m1.setPars({i*5+13:"0.,0.1"})
    m1.setPars({i*5+14:"0.,0.1"})
    Fit.perform()
    Fit.error("maximum 150 3.0 1-16")
    with open(outpara_name,'a')as f:
        for pa in range(1,17):
            try:
                print(pa,m1(pa).name,m1(pa).values[0],m1(pa).error[0],m1(pa).error[1],file=f)
            except:
                print(pa,m1(pa).name,m1(pa).values[0],-1,-1,file=f)
        print('*','CHIsqure',Fit.statistic,Fit.dof,file=f)
    m1.setPars({i*5+10:"0.,0.1"})
    m1.setPars({i*5+11:"0.,0.1"})
    Fit.perform()
    Fit.error("maximum 150 3.0 1-16")
    with open(outpara_name,'a')as f:
        for pa in range(1,17):
            try:
                print(pa,m1(pa).name,m1(pa).values[0],m1(pa).error[0],m1(pa).error[1],file=f)
            except:
                print(pa,m1(pa).name,m1(pa).values[0],-1,-1,file=f)
        print('*','CHIsqure',Fit.statistic,Fit.dof,file=f)
    for p in range(8):
        m1(p+7).frozen=True
    m1.setPars({15:"1.,0.1"})
    m1.setPars({16:"%s,-0.1"%crab_flux})
    Fit.perform()
    Fit.error("3.0 15")
    with open(outpara_name,'a')as f:
        for pa in range(1,17):
            try:
                print(pa,m1(pa).name,m1(pa).values[0],m1(pa).error[0],m1(pa).error[1],file=f)
            except:
                print(pa,m1(pa).name,m1(pa).values[0],-1,-1,file=f)
        print('*','CHIsqure',Fit.statistic,Fit.dof,file=f)
    #print>>f,' '
    #nmabs.append(m1(2).values[0])
    #psai.append([round(m1(7).values[0],3),round(2*m1(7).values[0]/(m1(7).error[1]-m1(7).error[0]),2)])
    #thta.append([round(m1(8).values[0],3),round(2*m1(8).values[0]/(m1(8).error[1]-m1(8).error[0]),2)])
    #phi.append([round(m1(9).values[0],3),round(2*m1(9).values[0]/(m1(9).error[1]-m1(9).error[0]),2)])
    m1.setPars({i*5+16:"%s,-1.0,0.,0.,1000.,1000,"%crab_flux})
    idx = np.array(d1.noticed)-1
    alfa = np.array(hxmtpsf.D_alfa)[idx]
    beta = np.array(hxmtpsf.D_beta)[idx]
    norm = np.array(d1.values)
    psf = np.array(m1.folded(1))
    Plot("data")
    yerr = Plot.yErr()
    yerr = np.array(yerr)
    yep = (np.abs(alfa)<8)&(np.abs(beta)<2.5)
    lists = np.c_[alfa[yep],beta[yep],norm[yep],yerr[yep],psf[yep]]
    np.save('%s_%s.npy'%(outnpy_name,lcnum),lists)
    AllModels.clear()
    AllData.clear()

'''
np.save('Roll_para_LE.npy',[psai,thta,phi])

psai,thta,phi = np.load('../ALL_old/Roll_para_LE.npy')

print psai
print thta
print phi

for lcnum in listt:
    AllData("%s%s.grp"%(outfile5,lcnum))
    AllData.ignore('81494-81700')
    AllData.ignore('75956-76200')
    AllData.ignore('50720-51300')
    AllData.ignore('64080-64120')
    AllModels+=temp
    NormLc=[]
    AlfaLc=[]
    erroLc=[]
    m1 = AllModels(1)
    d1 = AllData(1)
    i=0
    m1.setPars({1:"0 -1",2:"1. 0.1"})
    Fit.query = "yes"
    m1.setPars({i*5+3:" 1   -1"})
    m1.setPars({i*5+4:"%s    -0.001"%(rollt[lcnum])})
    m1.setPars({i*5+5:"83.633,-1.0"})
    m1.setPars({i*5+6:"22.015,-1.0"})
    m1.setPars({i*5+7:"%s,-1.0"%psai[lcnum][0]})
    m1.setPars({i*5+8:"%s,-1.0"%thta[lcnum][0]})
    m1.setPars({i*5+9:"%s,-1.0"%phi[lcnum][0]})
    m1.setPars({i*5+10:"260,1.0,0.,0.,1000.,1000,"})
    Fit.perform()
    idx = np.array(d1.noticed)-1
    alfa = np.array(hxmtpsf.D_alfa)[idx]
    beta = np.array(hxmtpsf.D_beta)[idx]
    m1.setPars({i*5+10:"260,-1.0,0.,0.,1000.,1000,"})
    norm = np.array(d1.values)
    psf = np.array(m1.folded(1))
    Plot("data")
    yerr = Plot.yErr()
    yerr = np.array(yerr)
    yep = (np.abs(alfa)<6.5)&(np.abs(beta)<2.)
    lists = np.c_[alfa[yep],beta[yep],norm[yep],yerr[yep],psf[yep]]
    np.save('Angle_normT_%s.npy'%(lcnum),lists)
    AllModels.clear()
    AllData.clear()


pa=np.load('Roll_para_LE.npy')
with open('para.txt','w')as f:
    for i in range(3):
        for j in range(3):
            print>>f,'%s\t%s'%(pa[i][j][0],pa[i][j][0]/pa[i][j][1])



for lcnum in listt:
    AllData("%s%s.grp"%(outfile5,lcnum))
    AllModels+=temp
    NormLc=[]
    AlfaLc=[]
    erroLc=[]
    m1 = AllModels(1)
    d1 = AllData(1)
    i=0
    m1.setPars({1:"0 -1",2:"%s 0.1"%10})
    Fit.query = "yes"
    m1.setPars({i*5+3:" 1   -1"})
    m1.setPars({i*5+4:"%s    -0.001"%(rollt[lcnum])})
    m1.setPars({i*5+5:"83.633,-1.0"})
    m1.setPars({i*5+6:"22.015,-1.0"})
    m1.setPars({i*5+7:"%s,-1.0"%0.})
    m1.setPars({i*5+8:"%s,-1.0"%0.})
    m1.setPars({i*5+9:"%s,-1.0"%0.})
    m1.setPars({i*5+10:"200,1.0,0.,0.,1000.,1000,"})
    Fit.perform()
    alfa = np.array(hxmtpsf.D_alfa)
    beta = np.array(hxmtpsf.D_beta)
    norm = np.array(d1.values)
    psf = np.array(m1.folded(1))
    Plot("data")
    yerr = Plot.yErr()
    yerr = np.array(yerr)
    yep = (np.abs(alfa)<6.2)&(np.abs(beta)<1.6)
    lists = np.c_[alfa[yep],beta[yep],norm[yep],yerr[yep],psf[yep]]
    np.save('Angle_normT_%s.npy'%(lcnum),lists)
    AllModels.clear()
    AllData.clear()
'''
##########################################################################################
########################gensrctxt
#gensrctxt(outfile4[:-16],outfile4[:-16])
##########################################################################################
logFile.close()
parFile.close()
time2 = time.time()
print('fit finished ',time2)
print("time expended: ",(time2-time1)/3600.0,' hours')

