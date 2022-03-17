#!/anaconda3/envs/FEALPy/bin python3.8
# -*- coding: utf-8 -*-
# ---
# @Software: PyCharm
# @File: psf_nurbs.py
# @Author: luoqi
# @Institution: Institute of High Energy Physics, CAS
# @E-mail: woshiluoqi123@outlook.com, luoqi@ihep.ac.cn
# @Site: 
# @Time: 3月 08, 2022
# ---
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

#-----------平面距离----------
def dis_surf(a0,b0,a1,b1):
    return ((a0-a1)**2 + (b0-b1)**2)**0.5

#-----------平方反比---------
def inv_sq(x):
    return 1/x**2
#---------------给一个点和一列点组，找到最近的四个dex和距离--------
def get_nearst4pts(x0,b0,data_x,data_y):
    dis_all = dis_surf(x0,b0,data_x,data_y)
    dex_array = np.argsort(dis_all)[0:4]
    dis_array = dis_all[dex_array]
    return dis_array,dex_array

def get_z_zerr(x0,b0,data_x,data_y,data_z,data_zerr):
    dis_array,dex_array = get_nearst4pts(x0,b0,data_x,data_y)
    nrst4z= data_z[dex_array]
    nrst4zerr = data_zerr[dex_array]
    if np.min(dis_array)==0:
        calc_z = data_z[dex_array[0]]
        calc_weight = data_zerr[dex_array[0]]
    else:
        calc_z = np.sum(inv_sq(nrst4zerr)*inv_sq(dis_array)*nrst4z)/np.sum(inv_sq(nrst4zerr)*inv_sq(dis_array))
        calc_weight = np.sum(inv_sq(nrst4zerr)*inv_sq(dis_array))/np.sum(inv_sq(dis_array))
    return calc_z,calc_weight

def nst_get_z_zerr(x0,y0,data_x,data_y,data_z,data_zerr):
    mindex = np.argmin((x0-data_x)**2+(y0-data_y)**2)
    calc_z = data_z[mindex]
    minzerr = data_zerr[mindex]
    calc_weight = inv_sq(minzerr)
    return calc_z,calc_weight

def read_data(instru,box_dex):
    data_x = np.load('%s_alpha_box%s.npy' % (instru,str(box_dex)))
    data_y = np.load('%s_beta_box%s.npy' % (instru, str(box_dex)))
    data_z = np.load('%s_cts_box%s.npy' % (instru, str(box_dex)))
    data_zerr = np.load('%s_cts_box%s.npy' % (instru, str(box_dex)))
    return data_x,data_y,data_z,data_zerr

def get_psf(x,y,data_x,data_y,data_z):
    # data_x = np.load('ME_alpha_box0.npy')
    # data_y = np.load('ME_beta_box0.npy')
    # data_z = np.load('ME_cts_box0.npy')
    data_zerr = np.ones_like(data_z)
    psf_value = np.zeros_like(x)
    if len(x)>1:
        for i in range(x.shape[0]):
            if len(x.shape)>1:
                for j in range(x.shape[1]):
                    calc_z, calc_weight = nst_get_z_zerr(x[i,j], y[i,j], data_x, data_y, data_z, data_zerr)
                    psf_value[i,j] = calc_z
            else:
                calc_z, calc_weight = nst_get_z_zerr(x[i], y[i], data_x, data_y, data_z, data_zerr)
                psf_value[i] = calc_z
    else:
        calc_z, calc_weight = nst_get_z_zerr(x, y, data_x, data_y, data_z, data_zerr)
        psf_value = calc_z
    # psf_value = psf_value + stand_psfmodel(x,y,instru=instru)
    return psf_value

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

def stand_psfmodel(delta_alfa0, delta_beta0,data_x,data_y,data_z, instru ='ME'):
    PI=3.14159265358979323846
    a_bound, b_bound = get_bound(instru)
    # roll_indx = box_dex#int(-parameter[1]/60+1)###用以确定机箱
    alfa_bound = np.abs(a_bound)#abound[roll_indx]
    beta_bound = np.abs(b_bound)#bbound[roll_indx]
    delta_alfa = np.abs(delta_alfa0)
    delta_beta = np.abs(delta_beta0)
    condi = np.logical_and(delta_alfa<alfa_bound,delta_beta<beta_bound)
    print(condi.shape)
    factor=np.where(condi,(1.0-np.tan(delta_alfa/180.0*PI)/np.tan(alfa_bound/180.0*PI))*\
        (1.0-np.tan(delta_beta/180.0*PI)/np.tan(beta_bound/180.0*PI)),np.zeros_like(delta_beta0))
    factor = factor / (np.sqrt(
        np.tan(delta_alfa / 180.0 * PI) * np.tan(delta_alfa / 180.0 * PI) + np.tan(delta_beta / 180.0 * PI) * np.tan(
            delta_beta / 180 * PI) + 1))
    return factor

def plot_psf(func, funcname,data_x,data_y,data_z,instru='ME'):
    fig = plt.figure(figsize=plt.figaspect(0.5))
    ax = fig.add_subplot(1, 1, 1)  # , projection='3d')
    a_bound,b_bound = get_bound(instru)
    alpha = np.arange(-a_bound, a_bound, 0.01)
    beta = np.arange(-b_bound, b_bound, 0.01)
    alpha, beta = np.meshgrid(alpha, beta)
    Z0 = func(alpha, beta,data_x,data_y,data_z)
    # Z1 = correct_psfmodel(alpha,beta)
    cmap = cm.YlOrRd
    ###cmap = cm.autumn
    ###cmap = cm.coolwarm
    levs = np.linspace(0.01, 1, 20)
    surf = ax.contourf(alpha, beta, Z0, levs,cmap=cmap)  # , rstride=1, cstride=1, cmap=cm.coolwarm,
    # linewidth=0, antialiased=False)
    plt.colorbar(surf, orientation='horizontal')
    fig.savefig('{:s}_{:s}_psf.png'.format(instru, funcname))
    return Z0

# plot_psf(stand_psfmodel, 'stand')
# correct_psf = plot_psf(get_psf, 'correct')


def plot_delta(data_x,data_y,data_z,instru='ME'):
    fig = plt.figure(figsize=plt.figaspect(0.5))
    ax = fig.add_subplot(1, 1, 1)  # , projection='3d')
    a_bound, b_bound = get_bound(instru)
    alpha = np.arange(-a_bound, a_bound, 0.01)
    beta = np.arange(-b_bound, b_bound, 0.01)
    alpha, beta = np.meshgrid(alpha, beta)
    Z0 = stand_psfmodel(alpha, beta,data_x,data_y,data_z)
    Z1 = nurbs_psf_module(alpha, beta, data_x, data_y, data_z)
    cmap = cm.YlOrRd
    ###cmap = cm.autumn
    ###cmap = cm.coolwarm
    levs = np.linspace(np.min(Z0-Z1), np.max(Z0-Z1), 20)
    surf = ax.contourf(alpha, beta, Z0 - Z1,levs, cmap=cmap)  # , rstride=1, cstride=1, cmap=cm.coolwarm,
    # linewidth=0, antialiased=False)
    plt.colorbar(surf, orientation='horizontal')
    fig.savefig('{:s}_psf_delta.png'.format(instru))

def nurbs_psf_module(alpha, beta, data_x, data_y, data_z):
    psf = (np.sqrt(2 * get_psf(alpha, beta, data_x, data_y, data_z) + 1) - 1) + stand_psfmodel(alpha, beta, data_x,
                                                                                         data_y, data_z)
    return psf

def plot_delta_minusstand(data_x,data_y,data_z,instru='ME'):
    fig = plt.figure(figsize=plt.figaspect(0.5))
    ax = fig.add_subplot(1, 1, 1)  # , projection='3d')
    a_bound, b_bound = get_bound(instru)
    alpha = np.arange(-a_bound, a_bound, 0.001)
    beta = np.arange(-b_bound, b_bound, 0.001)
    alpha, beta = np.meshgrid(alpha, beta)
    Z0 = stand_psfmodel(alpha, beta,data_x,data_y,data_z)
    Z1 = nurbs_psf_module(alpha, beta, data_x, data_y, data_z)
    # Z1 = np.exp(get_psf(alpha, beta,data_x,data_y,data_z))*stand_psfmodel(alpha, beta,data_x,data_y,data_z)
    cmap = cm.coolwarm
    surf = ax.contourf(alpha, beta, Z0 - Z1, cmap=cmap)  # , rstride=1, cstride=1, cmap=cm.coolwarm,
    # linewidth=0, antialiased=False)
    plt.colorbar(surf, orientation='horizontal')
    fig.savefig('{:s}_psf_delta.png'.format(instru))

def plot_psf_minustand(func, funcname,data_x,data_y,data_z,instru='ME'):
    fig = plt.figure(figsize=plt.figaspect(0.5))
    ax = fig.add_subplot(1, 1, 1)  # , projection='3d')
    a_bound,b_bound = get_bound(instru)
    alpha = np.arange(-a_bound, a_bound, 0.001)
    beta = np.arange(-b_bound, b_bound, 0.001)
    alpha, beta = np.meshgrid(alpha, beta)
    Z0 = func(alpha, beta,data_x,data_y,data_z)
    #Z0 = (np.sqrt(2 * get_psf(alpha, beta, data_x, data_y, data_z) + 1) - 1) + stand_psfmodel(alpha, beta, data_x,
    #                                                                                     data_y, data_z)
    #Z0 = (func(alpha, beta,data_x,data_y,data_z) * stand_psfmodel(alpha, beta,data_x,data_y,data_z))+stand_psfmodel(alpha, beta,data_x,data_y,data_z)
    # Z1 = correct_psfmodel(alpha,beta)
    cmap = cm.coolwarm
    surf = ax.contourf(alpha, beta, Z0, cmap=cmap)  # , rstride=1, cstride=1, cmap=cm.coolwarm,
    # linewidth=0, antialiased=False)
    plt.colorbar(surf, orientation='horizontal')
    fig.savefig('{:s}_{:s}_psf.png'.format(instru, funcname))
    return Z0

instru = 'ME'
psf_load = np.load('%s_psf_box0.npy'%(instru))
data_x = psf_load[:,0]
data_y = psf_load[:,1]
data_z = psf_load[:,2]
plot_psf(stand_psfmodel, 'stand',data_x,data_y,data_z)
# correct_psf = plot_psf(get_psf, 'correct',data_x,data_y,data_z)
# plot_delta(data_x,data_y,data_z,'ME')

correct_psf = plot_psf(nurbs_psf_module, 'correct',data_x,data_y,data_z)
plot_delta(data_x,data_y,data_z,'ME')