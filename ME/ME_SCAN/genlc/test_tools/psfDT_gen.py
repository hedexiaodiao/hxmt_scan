#!/hxmt/home/lick/soft/anaconda2/bin/python
'''
Copyright 2016, HXMT Ground System
All right reserved
File Name: hxmtpsf.py for sample

.......

Author: lick
Version: 0.1
Date:  2016-12-11
Email: lick@ihep.ac.cn
......
''' 
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
############read ra dec for scan stripe################
if len(sys.argv)<2:
    print "Need the config file!"
    cfg = "config.xml"
    #sys.exit(1)
else:
    cfg = sys.argv[1]
readcfg = loadDom(cfg)
#f = open("degree.dat","a")
infilestr = readcfg.getTagText("infilenamelist")
inpathstr = readcfg.getTagText("inpath")
instr = readcfg.getTagText("Instrument")
#print inpathstr,infilestr
inpathstr = inpathstr.strip()
evtfilestr = infilestr.split()[0]
infilestr2 = infilestr.split()[-1]
infile = (inpathstr+infilestr2).encode()
evtfile = (inpathstr+evtfilestr).encode()
print ("the scanning pointing file : ",infile)
#time.sleep(5) 
hdulist = pf.open(infile)
evtshdu = pf.open(evtfile)
t1 = hdulist[3].data
tb2 = evtshdu[1].data
qtime = t1.field(0)
etime = tb2.field(0)
tstart = etime[0]
tstop = etime[-1]
#q1_list = t1.field(1)[(qtime<=tstop) & (qtime>=tstart)]
#q2_list = t1.field(2)[(qtime<=tstop) & (qtime>=tstart)]
#q3_list = t1.field(3)[(qtime<=tstop) & (qtime>=tstart)]

q1_list = t1.field(1)[np.in1d(qtime,etime)]
q2_list = t1.field(2)[np.in1d(qtime,etime)]
q3_list = t1.field(3)[np.in1d(qtime,etime)]

instr = instr.strip()
instr = instr.split()[0]
instr = instr.encode()
#print instr,type(instr),len(instr),instr[0],instr[1]
if instr == "HE":
    long_bound=5.7 #means long_bound
    short_bound=1.1 #means short_bound
elif instr == "ME":
    long_bound=4.0 #means long_bound
    short_bound=1.0 #means short_bound
elif instr == "LE":
    long_bound=6.0 #means long_bound
    short_bound=1.6 #means short_bound
else:
    print "The instrment in config file is wrong."
    sys.exit(1)
print "The instrment and the bound :",instr, long_bound, short_bound
time.sleep(5)
mtx_list=[]

for i in range(0,q1_list.shape[0]):
    quat1 = [q1_list[i],q2_list[i],q3_list[i]]
    mtx_list.append(testroll.quat2mtx(quat1))

D_ra = []
D_dec= []
D_ra0= []
D_dec0= []


def hxmtpsf(engs, parameter, flux):
    PI=3.14159265358979323846
    flag=parameter[0]
    roll=parameter[1]/180.0*PI
    r_degree=parameter[1]
    ra=parameter[2]
    dec=parameter[3]
    idx = int(-(parameter[1]/60) + 1 )
    psai = twst['psai'][idx]
    theta=twst['theta'][idx]
    phi=twst['phi'][idx]
    NS=twst['NS'][idx]
    NL=twst['NL'][idx]
    rot_quat =  quat.Quat((psai,theta,phi))
    A = rot_quat.transform.T
    s_step=0
    Nflux = len(engs)-1
    x,y,z = testroll.sph2cart(ra*pi/180,dec*pi/180,1)
    xyz = np.array([x,y,z])
    dcm_b_f = np.array([[0,-sin(roll),cos(roll)],[0,-cos(roll),-sin(roll)],[1,0,0]]) 
    dcm_b_f0 = np.array([[0,-sin(np.pi/2.),cos(np.pi/2.)],[0,-cos(np.pi/2.),-sin(np.pi/2.)],[1,0,0]])
    ra_tplt=[]
    dec_tplt=[]
    dec0_tplt=[]
    ra0_tplt=[]

    for s_step in range(0,Nflux): 

        mtx = mtx_list[s_step]
        delta_long,delta_short,xr,yr,zr = testroll.quat2delta2_rot(mtx,dcm_b_f,xyz,A)# A?
        delta_long0,delta_short0,xr0,yr0,zr0 = testroll.quat2delta2_rot(mtx,dcm_b_f0,xyz,A)
        ra0_tplt.append(delta_long0)
        ra_tplt.append(delta_long)
        dec0_tplt.append(delta_short0)
        dec_tplt.append(delta_short)
        delta_long=fabs(delta_long)
        delta_short=fabs(delta_short)
        if delta_long<long_bound and delta_short<short_bound: 
            factor=(1.0-tan(delta_long/180.0*PI)/tan(long_bound/180.0*PI))*\
            (1.0-tan(delta_short/180.0*PI)/tan(short_bound/180.0*PI))
        else:
            factor=0
        factor=factor/(sqrt(tan(delta_long/180.0*PI)*tan(delta_long/180.0*PI)+tan(delta_short/180.0*PI)*tan(delta_short/180*PI)+1))
        if z !=0:
            factor=factor*zr*(1.0-NL*(xr/zr)*(xr/zr))*(1.0-NS*(yr/zr)*(yr/zr))

        flux[s_step]=factor
    D_ra[:]=ra_tplt
    D_ra0[:]=ra0_tplt
    D_dec[:]= dec_tplt
    D_dec0[:]= dec0_tplt

myModelParInfo=(
    "flag       \"\"      1     -10   -10    10   10   1", 
    "roll       deg    0     -90  -90   90  90  0.1",
    "ra         deg    0     0  0   360  360  0.5",
    "dec        deg    0     -90   -90    90   90   0.5")

