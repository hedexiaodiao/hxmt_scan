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
import psf_nurbs_v3

instru = "ME"
npy_name = "Angle"
lcnums = ['normC_0','normC_1','normC_2']



for i in range(3):
    psf_load = np.load('%s_psf_box%s.npy' % (instru, str(i)))
    data_x = psf_load[:, 0]
    data_y = psf_load[:, 1]
    data_z = psf_load[:, 2]

    # alpha = np.load('%s_alpha_box%s.npy' % (instru, str(i)))
    # beta = np.save('%s_beta_box%s' % (instru, str(i)))
    # norm = np.save('%s_cts_box%s' % (instru, str(i)))
    # yerr = np.save('%s_err_box%s' % (instru, str(i)))
    lists = np.transpose(np.load('%s_%s.npy'%(npy_name,lcnums[i])))
    #lists = np.c_[alfa[yep],beta[yep],norm[yep],yerr[yep],psf[yep]]
    print(lists.shape)
    alpha = lists[0]
    beta = lists[1]
    norm = lists[2]
    yerr = lists[3]
    psf = psf_nurbs_v3.nurbs_psf_module(alpha, beta, data_x, data_y, data_z)
    psf = psf*140#np.max(norm)
    #psf = lists[4]
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
    fig.savefig('fit_nurbs_%s_%s.png' % (npy_name,lcnums[i]))