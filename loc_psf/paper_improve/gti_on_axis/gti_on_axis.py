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

'''
To calibrate PSF with Crab in FOV alpha<0.8, beta<0.8.
No need to distinguish 3 box.
#  Crab, ra = 83.633, dec = 22.014
    H 1730-333, ra = 263.353, dec = -33.3888
    GX 354-0, ra = 262.989, dec = -33.835
1. calculate crt att with lc time
2. use Crab ra/dec to get alpha/beta list
3. return lc time dex
4. merge time for use 
'''


if len(sys.argv)<2:
    print("Need the config file!")
    cfg = "config.xml"
    #sys.exit(1)
else:
    cfg = sys.argv[1]

try:
    para_type = sys.argv[2]
except:
    para_type = 'all'

#---------------------load file-  -------------#
readcfg = loadDom(cfg)
infilestr = readcfg.getTagText("infilenamelist")
inpathstr = readcfg.getTagText("inpath")
instr = readcfg.getTagText("Instrument")
#print inpathstr,infilestr
inpathstr = inpathstr.strip()
evtfilestr = infilestr.split()[0]
infilestr2 = infilestr.split()[-1]
infile = (inpathstr+infilestr2)
evtfile = (inpathstr+evtfilestr)
print ("the scanning pointing file : ",infile)
instr = instr.strip()
instr = instr.split()[0]
instr = instr#.encode()

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
    print("The instrment in config file is wrong.")
    sys.exit(1)

#---------------------Crab ra/dec--------------#
ra =  83.633
dec = 22.014

alpha_lim = 0.8
beta_lim = 0.8

box_roll = [-60,0,60]
roll=box_roll[1]/180.0*np.pi
roll_indx = int(-box_roll[1]/60+1)

#--------------------read fits-----------------#
atthd = pf.open(infile)
lchd = pf.open(evtfile)
atdt = atthd[3].data
lcdt = lchd[1].data
qtime = atdt.field(0)
lctime = lcdt.field(0)
tstart = lctime[0]
tstop = lctime[-1]

#----------------quat with lc time------------#
q1_list = atdt.field(1)[np.in1d(qtime,lctime)]
q2_list = atdt.field(2)[np.in1d(qtime,lctime)]
q3_list = atdt.field(3)[np.in1d(qtime,lctime)]


#-------------calculate alpha/beta------------#
mtx_list = []
alpha_list = []
beta_list = []
for i in range(0,q1_list.shape[0]):
    quat1 = [q1_list[i],q2_list[i],q3_list[i]]
    mtx_list.append(testroll.quat2mtx(quat1))

psai = 0
theta = 0
phi = 0
rot_quat =  quat.Quat((psai,theta,phi))
x,y,z = testroll.sph2cart(ra*np.pi/180,dec*np.pi/180,1)
xyz = [x,y,z]
dcm_b_f_rot = rot_quat.transform.T


print('length: lighr curvefile, Quat file >>>',len(lctime),len(q1_list))

atthd.close()
lchd.close()
del atthd,lchd


mtx_list= []
accept_dex = []
accept_time = []
for i in range(0,q1_list.shape[0]):
    quat1 = [q1_list[i],q2_list[i],q3_list[i]]
    mtx = testroll.quat2mtx(quat1)
    mtx_list.append(mtx)
    delta_alfa0, delta_beta0, xr, yr, zr = testroll.quat2delta2_rot(mtx, dcm_b_f, xyz, dcm_b_f_rot)
    ###cosz = np.sqrt(zr ** 2 / (xr ** 2 + yr ** 2 + zr ** 2))
    delta_alfa = fabs(delta_alfa0 * cos(roll) - delta_beta0 * sin(roll))
    delta_beta = fabs(delta_alfa0 * sin(roll) + delta_beta0 * cos(roll))
    if delta_alfa<0.8 and delta_beta<0.8:
        accept_dex.append(i)
        accept_time.append(lctime[i])

merge_scal = 5
tem_accept_time = np.append([0],accept_time[1:])
gti_condi = np.greater(accept_time-tem_accept_time, 5)
gti_time_up = accept_time[gti_condi]
if len(gti_time_up)>1:
    gti_time_down = np.append(gti_time_up[1:],accept_time[-1])
elif len(gti_time_up)==1:
    gti_time_down = accept_time[-1]
else:
    print("No accept god time interval!!")

print(gti_time_up,'\n',gti_time_down)
###todo:直接保留time、att、lc文件写入fits便可，无须整理gti

