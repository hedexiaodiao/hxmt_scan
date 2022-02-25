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
from astropy.io import fits as pf
import sys
from plot_cald_psf_sl import hxmtpsf
from readxml import *
print("import success")

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

readcfg = loadDom(cfg)
infilestr = readcfg.getTagText("infilenamelist")
inpathstr = readcfg.getTagText("inpath")
instr = readcfg.getTagText("Instrument")
#print inpathstr,infilestr
inpathstr = inpathstr.strip()
evtfilestr_list = [infilestr.split()[0],infilestr.split()[1],infilestr.split()[2]]
infilestr2 = infilestr.split()[-1]
infile = (inpathstr+infilestr2)#.encode('UTF-8')
#print(infile)
#infile = (inpathstr+infilestr2).encode('GBK')
#print(infile)
#infile = (inpathstr+infilestr2)
#print(infile)
#.encode('UTF-8')
atthd = pf.open(infile)
atdt = atthd[3].data
qtime = atdt.field(0)

roll_list = [60.0,0.0,-60.0]
load_para = np.loadtxt(sys.argv[2])[:,0]
crab_flux = 282

for i in range(3):
    evtfile = (inpathstr + evtfilestr_list[i])
    print(evtfile)
    lchd = pf.open(evtfile)
    lcdt = lchd[1].data
    lctime = lcdt['Time']
    lc_flux = lcdt['counts']
    yerr = lcdt['error']
    tstart = lctime[0]
    tstop = lctime[-1]
    engs = lctime##for calc len only
    para_a = np.array([0.0,roll_list[i],83.633,22.015])
    ###print(para_a)
    ###print(load_para)
    para_b = load_para[int(i):int(9*(i+1))]
    ###print(para_b)
    parameter = np.append(para_a,para_b)#[flag,roll,ra,dec,psi,theta,phi,NS,NL,A,B,C,D]
    print(parameter)
    factor = hxmtpsf(engs, parameter, lc_flux)
    psf_flux = crab_flux*factor
    x_values = lctime
    psf = psf_flux
    norm = lc_flux
    sigma = (norm-psf)/yerr
    color_tem = 'black'
    fig = plt.figure(figsize=(30, 10))  # 10,20
    gs = gridspec.GridSpec(6, 3)
    ax1 = plt.subplot(gs[0:3, :])
    ax2 = plt.subplot(gs[3:5, :], sharex=ax1)
    ax3 = plt.subplot(gs[5:, :], sharex=ax1)
    gs.update(hspace=0)
    alpha = 0*norm###for convenient
    beta = 0*norm
    ax2.plot(x_values, alpha,'.', label='alpha')
    ax2.plot(x_values, beta,'.', label='beta')
    ax1.errorbar(x_values, norm, yerr, label='lc',fmt='_',ecolor=color_tem,mfc=color_tem, mfcalt=color_tem,mec=color_tem
                 ,markersize=6,elinewidth=1,capsize=0)
    ax1.plot(x_values, psf,'.', label='psf')
    ax3.errorbar(x_values, sigma, 1,label='sigma',fmt='_',ecolor=color_tem,mfc=color_tem, mfcalt=color_tem,mec=color_tem
                 ,markersize=6,elinewidth=1,capsize=0)
    ax1.legend()
    ax2.legend()
    ax3.legend()
    fig.savefig('%s_%s.png' % ('psf_box',i))