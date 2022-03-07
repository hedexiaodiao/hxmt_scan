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
    print('Need the config file!')
    cfg = 'config.xml'
else:
    cfg = sys.argv[1]
readcfg = loadDom(cfg)
infilestr = readcfg.getTagText('infilenamelist')
inpathstr = readcfg.getTagText('inpath')
instr = readcfg.getTagText('Instrument')
inpathstr = inpathstr.strip()
evtfile_list = [infilestr.split()[0],infilestr.split()[1],infilestr.split()[2]]
print('inpathstr:',inpathstr,'evt_list:',evtfile_list)
infilestr2 = infilestr.split()[-1]
infile = (inpathstr + infilestr2)###.encode()
###.encode()
print('the scanning pointing file : ', infile)

instr = instr.strip()
instr = instr.split()[0]
instr = instr###.encode()
if instr == 'HE':
    Talfa_bound = 5.7
    Tbeta_bound = 1.1
    dcm_b_f = np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]])
elif instr == 'ME':
    Talfa_bound = 4
    Tbeta_bound = 1
    dcm_b_f = np.array([[0, 0, -1], [0, 1, 0], [1, 0, 0]])
elif instr == 'LE':
    Talfa_bound = 6
    Tbeta_bound = 1.6
    dcm_b_f = np.array([[0, 0, 1], [0, 1, 0], [-1, 0, 0]])
else:
    print('The instrment in config file is wrong.')
    sys.exit(1)
print('The instrment and the bound :', instr, Talfa_bound, Tbeta_bound)
time.sleep(5)


D_alfa = []
D_beta = []

def hxmtpsf(engs, parameter, flux):
    PI = 3.14159
    flag = parameter[0]
    roll = (parameter[1] / 180) * PI
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

    #----------------------------gen mtx_list----------------
    box_index = int(-parameter[1] / 60 + 1)
    evtfile = (inpathstr + evtfile_list[box_index])
    print('parameter[0]:',parameter[0])
    print('parameter:',parameter)
    atthd = pf.open(infile)
    lchd = pf.open(evtfile)
    atdt = atthd[3].data
    lcdt = lchd[1].data
    qtime = atdt.field(0)
    lctime = lcdt.field(0)
    tstart = lctime[0]
    tstop = lctime[-1]
    q1_list = atdt.field(1)[np.in1d(qtime, lctime)]
    q2_list = atdt.field(2)[np.in1d(qtime, lctime)]
    q3_list = atdt.field(3)[np.in1d(qtime, lctime)]
    print('>>>>>>>>>>>>>>>>>>>>>>>>>>---->>>>>>>>>>>>>>>>>>>>>>>>>>')
    print(len(lctime), len(q1_list))
    atthd.close()
    lchd.close()
    del atthd
    del lchd
    mtx_list = []
    for i in range(0, q1_list.shape[0]):
        quat1 = [q1_list[i], q2_list[i], q3_list[i]]
        mtx_list.append(testroll.quat2mtx(quat1))
    print('Nflux,len_mtx:',Nflux,len(mtx_list))
    #---------------------------------------------------

    for s_step in range(0, Nflux):
        mtx = mtx_list[s_step]
        delta_alfa0, delta_beta0, xr, yr,zr = testroll.quat2delta2_rot(mtx, dcm_b_f, xyz, dcm_b_f_rot)
        cosz = np.sqrt(zr ** 2 / (xr ** 2 + yr ** 2 + zr ** 2))
        delta_alfa = fabs(delta_alfa0 * cos(roll) - delta_beta0 * sin(roll))
        delta_beta = fabs(delta_alfa0 * sin(roll) + delta_beta0 * cos(roll))
        alfa_tplt.append(delta_alfa0 * cos(roll) - delta_beta0 * sin(roll))
        beta_tplt.append(delta_alfa0 * sin(roll) + delta_beta0 * cos(roll))
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
