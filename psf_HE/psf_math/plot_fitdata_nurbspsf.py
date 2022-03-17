#!/anaconda3/envs/FEALPy/bin python3.8
# -*- coding: utf-8 -*-
# ---
# @Software: PyCharm
# @File: plot_fitdata_nurbspsf.py
# @Author: luoqi
# @Institution: Institute of High Energy Physics, CAS
# @E-mail: woshiluoqi123@outlook.com, luoqi@ihep.ac.cn
# @Site: 
# @Time: 3月 17, 2022
# ---

import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec

npy_name = "Angle"
lcnums = ['normC_0','normC_1','normC_2']

for i in range(3):
    lists = np.transpose(np.load('%s_%s.npy'%(npy_name,lcnums[i])))
#lists = np.c_[alfa[yep],beta[yep],norm[yep],yerr[yep],psf[yep]]
    print(lists.shape)
    alpha = lists[0]
    beta = lists[1]
    norm = lists[2]
    yerr = lists[3]
    psf = lists[4]
    sigma = (norm-psf)/yerr
    print(len(alpha))
    x_values = np.arange(0, len(alpha))
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
    fig.savefig('%s_%s.png' % (npy_name,lcnums[i]))