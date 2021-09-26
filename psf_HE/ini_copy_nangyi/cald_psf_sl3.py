from xspec import *
from math import *
import numpy as np
from astropy.io import fits as pf
from readxml import *
import time
import sys
import Quat as quat
import testroll
if len(sys.argv) < 2:
    print 'Need the config file!'
    cfg = 'config.xml'
else:
    cfg = sys.argv[1]
readcfg = loadDom(cfg)
infilestr = readcfg.getTagText('infilenamelist')
inpathstr = readcfg.getTagText('inpath')
instr = readcfg.getTagText('Instrument')
inpathstr = inpathstr.strip()
evtfilestr = infilestr.split()
infilestr2 = infilestr.split()

qfile = [(inpathstr+i).encode() for i in infilestr2][3:]
evtfile = [(inpathstr+i).encode() for i in infilestr2][:3]
print ("the scanning pointing file : ",qfile)
mtx_list=[]
q1s=[]
q2s=[]
q3s=[]
for i in range(3):
    atthd = pf.open(qfile[i])
    lchd = pf.open(evtfile[i])
    atdt = atthd[3].data
    lcdt = lchd[1].data
    qtime = atdt.field(0)
    lctime = lcdt.field(0)

    q1_list = atdt.field(1)[np.in1d(qtime,lctime)]
    q2_list = atdt.field(2)[np.in1d(qtime,lctime)]
    q3_list = atdt.field(3)[np.in1d(qtime,lctime)]
    q1s.append(q1_list);q2s.append(q2_list);q3s.append(q3_list)
    mtxs=[]
    for i in range(0,q1_list.shape[0]):
        quat1 = [q1_list[i],q2_list[i],q3_list[i]]
        mtxs.append(testroll.quat2mtx(quat1))

    atthd.close()
    lchd.close()
    del atthd,lchd
    mtx_list.append(mtxs)

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


D_alfa = []
D_beta = []

def hxmtpsf(engs, parameter, flux):
    PI = 3.14159
    flag = parameter[0]
    roll = (parameter[1] / 180) * PI
    roll_index=int(-parameter[1]/60+1)
    r_degree = parameter[1]
    ra = parameter[2]
    dec = parameter[3]
    psai = parameter[4]
    theta = parameter[5]
    phi = parameter[6]
    alfa_bound = Talfa_bound + parameter[7]
    beta_bound = Tbeta_bound + parameter[8]
    (pa, pb, pc, pd) = (parameter[9], parameter[10], parameter[11], parameter[12])
    s_step = 0
    Nflux = len(engs) - 1
    rot_quat = quat.Quat((psai, theta, phi))
    (x, y, z) = testroll.sph2cart(ra * pi / 180, dec * pi / 180, 1)
    xyz = [x, y, z]
    dcm_b_f_rot = rot_quat.transform.T
    alfa_tplt = []
    beta_tplt = []
    #print Nflux,roll_index,len(mtx_list[roll_index])
    for s_step in range(0, Nflux):
        mtx = mtx_list[roll_index][s_step]
        #delta_alfa0, delta_beta0, xr, yr,zr = testroll.quat2delta2_rot(mtx, dcm_b_f, xyz, dcm_b_f_rot)
        xr,yr,zr=np.dot([[cos(roll),-sin(roll),0],[sin(roll),cos(roll),0],[0,0,1.]],np.dot(dcm_b_f,np.dot(mtx,xyz)))
        delta_alfa0, delta_beta0 = atan2(xr,zr),atan2(yr,zr)
        delta_alfa=fabs(delta_alfa0)
        delta_beta=fabs(delta_beta0)
        alfa_tplt.append(delta_alfa0)
        beta_tplt.append(delta_beta0)
        if delta_alfa < alfa_bound and delta_beta < beta_bound:
            factor = (1 - tan((delta_alfa / 180) * PI) / tan((alfa_bound / 180) * PI)) * (1 - tan((delta_beta / 180) * PI) / tan((beta_bound / 180) * PI))
        else:
            factor = 0
        factor = factor / sqrt(tan((delta_alfa / 180) * PI) * tan((delta_alfa / 180) * PI) + tan((delta_beta / 180) * PI) * tan((delta_beta / 180) * PI) + 1)
        rt = pa * delta_alfa ** 2 * delta_beta ** 2 + pb * delta_alfa ** 2 + pc * delta_beta ** 2 + pd
        flux[s_step] = factor * rt

    D_alfa[:] = alfa_tplt
    D_beta[:] = beta_tplt

myModelParInfo = (
    'flag       ""      1     -10   -10    10   10   1',
    'roll       deg    0     -90  -90   90  90  0.1',
    'ra         deg    0     0  0   360  360  0.5',
    'dec        deg    0     -90   -90    90   90   0.5',
    'psai         deg    0     -25  -25   25  25  -0.001',
    'theta         deg    0     -25  -25   25  25  -0.001',
    'phi         deg    0     -25  -25   25  25  -0.001',
    'NS         deg    0     -25  -25   25  25  -0.001',
    'NL         deg    0     -25  -25   25  25  -0.001',
    'A         deg    0     -100000  -100000   100000  100000  -0.001',
    'B         deg    0     -100000  -100000   100000  100000  -0.001',
    'C         deg    0     -100000  -100000   100000  100000  -0.001',
    'D         deg    1.     -100000  -100000   100000  100000  -0.001')
