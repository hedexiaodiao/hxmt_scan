#!/anaconda3/envs/FEALPy/bin python3.8
# -*- coding: utf-8 -*-
# ---
# @Software: PyCharm
# @File: plot_fitdata.py
# @Author: luoqi
# @Institution: Institute of High Energy Physics, CAS
# @E-mail: woshiluoqi123@outlook.com, luoqi@ihep.ac.cn
# @Site: 
# @Time: 1月 21, 2022
# ---
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec

npy_name = "A02"
lcnums = ['me_b0','me_b1','me_b2']

import matplotlib
matplotlib.use("agg")
from math import *
import numpy as np
from astropy.io import fits as pf
from readxml import *
import time
import sys
import Quat as quat
#from testroll import *
import testroll
from xml.etree.ElementTree import ElementTree,Element
import matplotlib.pyplot as plt

alpha_lim = 1.6
alpha_uplim = [alpha_lim,alpha_lim,alpha_lim]
alpha_downlim = [-alpha_lim,-alpha_lim,-alpha_lim]
beta_lim = 0.8
beta_uplim = [beta_lim,beta_lim,beta_lim]
beta_downlim = [-beta_lim,-beta_lim,-beta_lim]


def stand_model(data_x, data_y, a_bound, b_bound):
    data_x = np.array(data_x)
    data_y = np.array(data_y)
    # ----------standard model---------
    print(data_x.shape,data_y.shape,data_x,data_y)
    condi = np.logical_and(data_x < a_bound, data_y < b_bound)
    rad_data_x = np.deg2rad(np.abs(data_x))
    rad_data_y = np.deg2rad(np.abs(data_y))
    rad_a_bound = np.deg2rad(a_bound)
    rad_b_bound = np.deg2rad(b_bound)
    factor = np.where(condi, (1.0 - np.tan(rad_data_x) / np.tan(rad_a_bound)) * \
                      (1.0 - np.tan(rad_data_y) / np.tan(rad_b_bound)),
                      np.zeros_like(rad_data_x))
    factor = factor / (np.sqrt(np.tan(rad_data_x) ** 2 + np.tan(rad_data_y) ** 2 + 1))
    return factor


def get_bound(instru):
    if instru == 'HE':
        Talfa_bound = 5.7
        Tbeta_bound = 1.1
        dcm_b_f = np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]])
    elif instru == 'ME':
        Talfa_bound = 4
        Tbeta_bound = 1
        dcm_b_f = np.array([[0, 0, -1], [0, 1, 0], [1, 0, 0]])
    elif instru == 'LE':
        Talfa_bound = 6
        Tbeta_bound = 1.6
        dcm_b_f = np.array([[0, 0, 1], [0, 1, 0], [-1, 0, 0]])
    return Talfa_bound, Tbeta_bound


def plot_origindata(norm, yerr, x_values, alpha, beta, instru):
    a_bound, b_bound = get_bound(instru)
    factor = stand_model(alpha, beta, a_bound, b_bound)
    psf = factor * np.max(cts) / np.max(factor)

    sigma = (cts - psf) / yerr
    print(len(alpha))
    ###x_values = np.arange(0, len(alpha))
    print(x_values)
    '''
    fig = plt.figure(figsize=plt.figaspect(0.5))
    ax = fig.add_subplot(1, 1, 1)  # , projection='3d')
    plt.plot(x_values, alpha,label='alpha')
    plt.plot(x_values, beta,label='beta')
    plt.errorbar(x_values, norm, yerr,label='lc')
    plt.plot(x_values,psf,label='psf')
    plt.legend()
    '''
    color_tem = 'black'
    fig = plt.figure(figsize=(30, 10))  # 10,20
    gs = gridspec.GridSpec(6, 3)
    ax1 = plt.subplot(gs[0:3, :])
    ax2 = plt.subplot(gs[3:5, :], sharex=ax1)
    ax3 = plt.subplot(gs[5:, :], sharex=ax1)
    gs.update(hspace=0)
    ax2.plot(x_values, alpha, '.', label='alpha')
    ax2.plot(x_values, beta, '.', label='beta')
    ax1.errorbar(x_values, norm, yerr, label='lc', fmt='_', ecolor=color_tem, mfc=color_tem, mfcalt=color_tem,
                 mec=color_tem
                 , markersize=6, elinewidth=1, capsize=0)
    ax1.plot(x_values, psf, '.', label='psf')
    ax3.errorbar(x_values, sigma, 1, label='sigma', fmt='_', ecolor=color_tem, mfc=color_tem, mfcalt=color_tem,
                 mec=color_tem
                 , markersize=6, elinewidth=1, capsize=0)
    ax1.legend()
    ax2.legend()
    ax3.legend()
    fig.savefig('%s_alim%s_blim%s.png' % (evtfile[:-5],alpha_lim,beta_lim))

if len(sys.argv)<2:
    print("Need the config file!")
    cfg = "config.xml"
    #sys.exit(1)
else:
    cfg = sys.argv[1]

# year = sys.argv[2]

if len(sys.argv)==4:
    select_box_dex = int(sys.argv[3])
else:
    select_box_dex = 1

#---------------make config------------------------#
def gen_config(att_list,lc_list,fits_dir,xmlpath='./config_me.xml',instru='LE'):
    exec_name = './command_'+instru.lower()+'_'+time.strftime("%y%m%d")+'.sh'

    tree = ElementTree()
    tree.parse(xmlpath)
    pathnodes = tree.findall("PATH/outpath")
    outnodes = tree.findall("Outfile_Name/outfilename")
    filenodes = tree.findall("Infile_Name/infilename/infilenamelist")
    instrunodes = tree.findall("Inst_Info/Instrument")
    outfile = fits_dir + "/config_%s.xml"%(instru.lower())
    print(outfile)
    filenodes[0].text = "\n  " + lc_list[0]+" \n"
    filenodes[1].text = "\n  " + lc_list[1]+" \n"
    filenodes[2].text = "\n  " + lc_list[2]+" \n"
    filenodes[3].text = '\n  ' + att_list[0] + ' \n'
    pathnodes[0].text = '\n  ' + fits_dir + '  \n'
    instrunodes[0].text = '\n ' +instru + ' \n'
    outnodes[0].text = '\n '+'/GAL_he_small \n'

    tree.write(outfile, encoding="utf-8", xml_declaration=True)
    with open(exec_name,'a+')as f:
        try:
            print('sh ./sub_task_psf.sh %s'%(outfile),file=f)
        except:
            print('cannot write down command!!!')

#---------------------load file-  -------------#
readcfg = loadDom(cfg)
infilestr = readcfg.getTagText("infilenamelist")
inpathstr = readcfg.getTagText("inpath")
outpathstr = readcfg.getTagText("outpath")
instr = readcfg.getTagText("Instrument")
#print inpathstr,infilestr
inpathstr = inpathstr.strip()
outpathstr = outpathstr.strip()
evtfilestr0 = infilestr.split()[0]
evtfilestr1 = infilestr.split()[1]
evtfilestr2 = infilestr.split()[2]
infilestr2 = infilestr.split()[-1]
infile = (inpathstr+infilestr2)
evtfile0 = (inpathstr+evtfilestr0)
evtfile1 = (inpathstr+evtfilestr1)
evtfile2 = (inpathstr+evtfilestr2)
evtfile_list = [evtfile0,evtfile1,evtfile2]
dir_front = outpathstr

evtfile = evtfile_list[select_box_dex]

print ("the scanning pointing file : ",infile)
print("the lc file :", evtfile_list)
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


box_roll = [60,0,-60]
roll=box_roll[select_box_dex]/180.0*np.pi
roll_indx = int(-box_roll[select_box_dex]/60+1)
lc_list = []

for select_box_dex in range(3):
    evtfile = evtfile_list[select_box_dex]
    # --------------------read fits-----------------#
    atthd = pf.open(infile)
    lchd = pf.open(evtfile)
    atdt = atthd[3].data
    lcdt = lchd[1].data
    qtime = atdt.field(0)
    lctime = lcdt.field(0)
    tstart = lctime[0]
    tstop = lctime[-1]

    # ----------------quat with lc time------------#
    q1_list = atdt.field(1)[np.in1d(qtime, lctime)]
    q2_list = atdt.field(2)[np.in1d(qtime, lctime)]
    q3_list = atdt.field(3)[np.in1d(qtime, lctime)]

    # -------------calculate alpha/beta------------#
    mtx_list = []
    alpha_list = []
    beta_list = []
    for i in range(0, q1_list.shape[0]):
        quat1 = [q1_list[i], q2_list[i], q3_list[i]]
        mtx_list.append(testroll.quat2mtx(quat1))

    psai = 0
    theta = 0
    phi = 0
    rot_quat = quat.Quat((psai, theta, phi))
    x, y, z = testroll.sph2cart(ra * np.pi / 180, dec * np.pi / 180, 1)
    xyz = [x, y, z]
    dcm_b_f_rot = rot_quat.transform.T

    print('length: lighr curvefile, Quat file >>>', len(lctime), len(q1_list))

    atthd.close()
    lchd.close()
    del atthd, lchd

    #------------------select--------------
    mtx_list = []
    accept_dex = []
    accept_time = []
    accept_alpha = []
    accept_beta = []
    for i in range(0,q1_list.shape[0]):
        quat1 = [q1_list[i],q2_list[i],q3_list[i]]
        mtx = testroll.quat2mtx(quat1)
        mtx_list.append(mtx)
        delta_alfa0, delta_beta0, xr, yr, zr = testroll.quat2delta2_rot(mtx, dcm_b_f, xyz, dcm_b_f_rot)
        ###cosz = np.sqrt(zr ** 2 / (xr ** 2 + yr ** 2 + zr ** 2))

        roll = box_roll[select_box_dex] / 180.0 * np.pi
        delta_alfa_box = delta_alfa0 * cos(roll) - delta_beta0 * sin(roll)
        delta_beta_box = delta_alfa0 * sin(roll) + delta_beta0 * cos(roll)
        if delta_alfa_box < alpha_uplim[select_box_dex] and delta_alfa_box > alpha_downlim[select_box_dex] and delta_beta_box < beta_uplim[select_box_dex] and delta_beta_box > beta_downlim[select_box_dex]:
            flag = 1
            #tot_flag = np.logical_and(np.logical_and(flag[0],flag[1]),flag[2])
            #if tot_flag:
            accept_dex.append(i)
            accept_time.append(lctime[i])
            accept_alpha.append(delta_alfa0)
            accept_beta.append(delta_beta0)

    boxfile = pf.open('%s' % (evtfile))
    t = boxfile[1].data.field(0)
    cts = boxfile[1].data.field(1)
    err = boxfile[1].data.field(2)
    ###bkg = boxfile[1].data.field('bkg')
    mask = np.in1d(t, accept_time)
    print("before select:",len(t))
    print("after select:",len(accept_time))
    t = t[mask]
    cts = cts[mask]
    err = err[mask]
    plot_origindata(norm=cts, yerr=err, x_values=t, alpha=accept_alpha, beta=accept_beta, instru=instr)


# att_list = [infile]
# gen_config(att_list,lc_list,dir_front)


# fig = plt.figure(figsize=plt.figaspect(0.5))
# ax = fig.add_subplot(1, 1, 1)  # , projection='3d')
# plt.plot(np.array(accept_alpha), np.array(accept_beta),"*")
# fig.savefig('%s_data_position.png'%(instr))

