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

infilestr = readcfg.getTagText("infilenamelist")
inpathstr = readcfg.getTagText("inpath")
instr = readcfg.getTagText("Instrument")
inpathstr = inpathstr.strip()
evtfilestr = infilestr.split()[0]
infilestr2 = infilestr.split()[-1]
infile = (inpathstr+infilestr2).encode()
evtfile = (inpathstr+evtfilestr).encode()
print ("the scanning pointing file : ",infile)
########################################################

Twst = {'psai':[0.097,-0.022,-0.061],'theta':[0.043,0.131,0.111],'phi':[-0.01,0.31,0.20],'NS':[0.,29.4,0.],'NL':[31.3,525.5,6.37]}

########################################################
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
    ra_bound=5.7 #means long_bound
    dec_bound=1.1 #means short_bound
elif instr == "ME":
    dec_bound=4.0 #means long_bound
    ra_bound=1.0 #means short_bound
elif instr == "LE":
    dec_bound=6.0 #means long_bound
    ra_bound=1.6 #means short_bound
else:
    print "The instrment in config file is wrong."
    sys.exit(1)
#print "The instrment and the bound :",instr, long_bound, short_bound
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
    psai = Twst['psai'][idx]
    theta=Twst['theta'][idx]
    phi=Twst['phi'][idx]
    NS=Twst['NS'][idx]
    NL=Twst['NL'][idx]
    s_step=0
    Nflux = len(engs)-1
    x,y,z = testroll.sph2cart(ra*pi/180,dec*pi/180,1)
    xyz = np.array([x,y,z])
    rot_quat =  quat.Quat((psai,theta,phi))
    dcm_b_f_rot = rot_quat.transform.T
    dcm_b_f = np.array([[0,1,0],[0,0,1],[1,0,0]]) 
    dcm_b_f0 = np.array([[0,-sin(np.pi/2.),cos(np.pi/2.)],[0,-cos(np.pi/2.),-sin(np.pi/2.)],[1,0,0]])
    ra_tplt=[]
    dec_tplt=[]
    dec0_tplt=[]
    ra0_tplt=[]

    for s_step in range(0,Nflux): 
        mtx = mtx_list[s_step]
        delta_ra0,delta_dec0,xr,yr,zr = testroll.quat2delta2_rot(mtx,dcm_b_f,xyz,dcm_b_f_rot)
        delta_ra=fabs(delta_ra0*cos(roll)-delta_dec0*sin(roll))
        delta_dec=fabs(delta_ra0*sin(roll)+delta_dec0*cos(roll))
        ra_tplt.append(delta_ra0*cos(roll)-delta_dec0*sin(roll))
        dec_tplt.append(delta_ra0*sin(roll)+delta_dec0*cos(roll))
        ra0_tplt.append(delta_ra0*cos(0.)-delta_dec0*sin(0.))
        dec0_tplt.append(delta_ra0*sin(0.)+delta_dec0*cos(0.))
        ra_crt0 = delta_ra0*cos(0.)-delta_dec0*sin(0.)
        dec_crt0 = delta_ra0*sin(0.)+delta_dec0*cos(0.)
        dec_crt= delta_ra0*sin(roll)+delta_dec0*cos(roll)
        ra_crt = delta_ra0*cos(roll)-delta_dec0*sin(roll)

        if delta_ra<ra_bound and delta_dec<dec_bound: 
            factor=(1.0-tan(delta_ra/180.0*PI)/tan(ra_bound/180.0*PI))*\
            (1.0-tan(delta_dec/180.0*PI)/tan(dec_bound/180.0*PI))
        else:
            factor=0
        factor=factor/(sqrt(tan(delta_ra/180.0*PI)*tan(delta_ra/180.0*PI)+tan(delta_dec/180.0*PI)*tan(delta_dec/180*PI)+1))
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

