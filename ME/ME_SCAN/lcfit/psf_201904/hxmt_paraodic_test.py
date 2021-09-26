from xspec import *
from math import *
import numpy as np
from astropy.io import fits as pf
from readxml import *
import time 
import sys
import Quat as quat
#from testroll import *
import testroll 

if len(sys.argv)<2:
    print "Need the config file!"
    cfg = '/sharefs/hbkg/user/saina/'+raw_input("config file: /sharefs/hbkg/user/saina/(1)*/(2)*/(3)*/ME/config_me_small.xml: ")+'/ME/config_me_small.xml'
    #sys.exit(1)
else:
    cfg = sys.argv[1]

try:
    para_type = sys.argv[2]
except:
    para_type = 'all'

readcfg = loadDom(cfg)
infilestr = readcfg.getTagText("infilenamelist")
inpathstr = readcfg.getTagText("inpath")
instr = readcfg.getTagText("Instrument")
#print inpathstr,infilestr
inpathstr = inpathstr.strip()###inpath
evtfilestr = infilestr.split()[0]
infilestr2 = infilestr.split()[-1]###Att file
infile = (infilestr2).encode()###
evtfile = (evtfilestr).encode()###(inpathstr+evtfilestr).encode()
print ("the scanning pointing file : ",infile)

atthd = pf.open(infile)
lchd = pf.open(evtfile)
atdt = atthd[3].data
lcdt = lchd[1].data
qtime = atdt.field(0)
lctime = lcdt.field(0)
tstart = lctime[0]
tstop = lctime[-1]
#q1_list = t1.field(1)[(qtime<=tstop) & (qtime>=tstart)]
#q2_list = t1.field(2)[(qtime<=tstop) & (qtime>=tstart)]
#q3_list = t1.field(3)[(qtime<=tstop) & (qtime>=tstart)]

q1_list = atdt.field(1)[np.in1d(qtime,lctime)]
q2_list = atdt.field(2)[np.in1d(qtime,lctime)]
q3_list = atdt.field(3)[np.in1d(qtime,lctime)]

print '>>>>>>>>>>>>>>>>>>>>>>>>>>---->>>>>>>>>>>>>>>>>>>>>>>>>>'
print len(lctime),len(q1_list)

atthd.close()
lchd.close()
del atthd,lchd

instr = instr.strip()
instr = instr.split()[0]
instr = instr.encode()
#print instr,type(instr),len(instr),instr[0],instr[1]
if instr == "HE":
    Talfa_bound=5.7 #means long_bound
    Tbeta_bound=1.1 #means short_bound
    dcm_b_f = np.array([[0,1,0],[0,0,1],[1,0,0]])
elif instr == "ME":
    Talfa_bound=4.0 #means long_bound
    Tbeta_bound=1.0 #means short_bound
    dcm_b_f = np.array([[0,0,-1],[0,1,0],[1,0,0]])
elif instr == "LE":
    Talfa_bound=6.0 #means long_bound
    Tbeta_bound=1.6 #means short_bound
    dcm_b_f = np.array([[0,0,1],[0,1,0],[-1,0,0]])
else:
    print "The instrment in config file is wrong."
    sys.exit(1)

print "The instrment and the bound :",instr, Talfa_bound, Tbeta_bound
time.sleep(5)
mtx_list=[]
for i in range(0,q1_list.shape[0]):
    quat1 = [q1_list[i],q2_list[i],q3_list[i]]
    mtx_list.append(testroll.quat2mtx(quat1))

D_alfa = []
D_beta= []
#paras = np.loadtxt('/hxmt/work/HXMT_scan_data/psfcrt/psf_201811/Parabolic_%s_18.txt'%instr)
paras = np.loadtxt('/sharefs/hbkg/data/SCAN/PSF/psf_201904/Parabolic_%s_%s_201904.txt'%(instr,para_type))###
#paras = np.loadtxt('/sharefs/hbkg/data/SCAN/PSF/psf_201904/Parabolic_LE_2018_201904.txt')

Tpsai=paras[::9,0]
Ttheta=paras[1::9,0]
Tphi=paras[2::9,0]
abound = Talfa_bound + paras[3::9,0]
bbound = Tbeta_bound + paras[4::9,0]
Tpa,Tpb,Tpc,Tpd = paras[5::9,0],paras[6::9,0],paras[7::9,0],paras[8::9,0]
print Tpsai,Tpc
def hxmtpsf(engs, parameter, flux):
    PI=3.14159265358979323846
    flag=parameter[0]
    roll=parameter[1]/180.0*PI
    roll_indx = int(-parameter[1]/60+1)
    r_degree=parameter[1]
    ra=parameter[2]
    dec=parameter[3]
    psai=Tpsai[roll_indx]
    theta=Ttheta[roll_indx]
    phi=Tphi[roll_indx]
    alfa_bound = abound[roll_indx]
    beta_bound = bbound[roll_indx]
    pa,pb,pc,pd = Tpa[roll_indx],Tpb[roll_indx],Tpc[roll_indx],Tpd[roll_indx]
    s_step=0
    Nflux = len(engs)-1
    rot_quat =  quat.Quat((psai,theta,phi))
    x,y,z = testroll.sph2cart(ra*pi/180,dec*pi/180,1)
    xyz = [x,y,z]
    dcm_b_f_rot = rot_quat.transform.T
    alfa_tplt=[]
    beta_tplt=[]
    for s_step in range(0,Nflux): 
        mtx = mtx_list[s_step]
        delta_alfa0,delta_beta0,xr,yr,zr = testroll.quat2delta2_rot(mtx,dcm_b_f,xyz,dcm_b_f_rot)
        cosz = np.sqrt(zr**2/(xr**2+yr**2+zr**2))
        delta_alfa=fabs(delta_alfa0*cos(roll)-delta_beta0*sin(roll))
        delta_beta=fabs(delta_alfa0*sin(roll)+delta_beta0*cos(roll))
        alfa_tplt.append(delta_alfa0*cos(roll)-delta_beta0*sin(roll))
        beta_tplt.append(delta_alfa0*sin(roll)+delta_beta0*cos(roll))
        if delta_alfa<alfa_bound and delta_beta<beta_bound: 
            factor=(1.0-tan(delta_alfa/180.0*PI)/tan(alfa_bound/180.0*PI))*\
            (1.0-tan(delta_beta/180.0*PI)/tan(beta_bound/180.0*PI))
        else:
            factor=0
        factor= factor/(sqrt(tan(delta_alfa/180.0*PI)*tan(delta_alfa/180.0*PI)+tan(delta_beta/180.0*PI)*tan(delta_beta/180*PI)+1))
        rt = pa * delta_alfa**2 * delta_beta**2 + pb*delta_alfa**2 +pc*delta_beta**2 + pd
        flux[s_step]=factor*rt
    D_alfa[:]=alfa_tplt
    D_beta[:]=beta_tplt

myModelParInfo=(
    "flag       \"\"      1     -10   -10    10   10   1", 
    "roll       deg    0     -90  -90   90  90  0.1",
    "ra         deg    0     0  0   360  360  0.5",
    "dec        deg    0     -90   -90    90   90   0.5") 

