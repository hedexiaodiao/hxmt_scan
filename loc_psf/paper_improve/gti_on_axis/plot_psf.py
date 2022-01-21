from xspec import *
from math import *
import numpy as np
from astropy.io import fits as pf
###from readxml import *
import time
import sys
import Quat as quat
#from testroll import *
import testroll
import matplotlib.pyplot as plt
from matplotlib import cm

instr = "LE"
#instr = instr.strip()
#instr = instr.split()[0]
#instr = instr.encode()
print("instr:",instr)
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
    print("The instrment in config file is wrong.")
    sys.exit(1)

print("The instrment and the bound :",instr, Talfa_bound, Tbeta_bound)



#paras = np.loadtxt('/hxmt/work/HXMT_scan_data/psfcrt/psf_201811/Parabolic_%s_18.txt'%instr)
###paras = np.loadtxt('./Parabolic_%s_%s_201904.txt'%(instr,str(2018)))
def psf(Astr):
    paras = np.loadtxt('./Parabolic_%s_%s_210825.txt'%(instr,str(Astr)))
    #paras = np.loadtxt('/sharefs/hbkg/data/SCAN/PSF/psf_201904/Parabolic_LE_2018_201904.txt')
    
    Tpsai=paras[::9,0]
    Ttheta=paras[1::9,0]
    Tphi=paras[2::9,0]
    abound = Talfa_bound + paras[3::9,0]
    bbound = Tbeta_bound + paras[4::9,0]
    Tpa,Tpb,Tpc,Tpd = paras[5::9,0],paras[6::9,0],paras[7::9,0],paras[8::9,0]
    print(Tpsai,Tpc)
    
    def stand_psfmodel(delta_alfa0, delta_beta0, box_dex=0):
        PI=3.14159265358979323846
        roll_indx = box_dex#int(-parameter[1]/60+1)###用以确定机箱
        alfa_bound = np.abs(Talfa_bound + 0)#abound[roll_indx]
        beta_bound = np.abs(Tbeta_bound + 0)#bbound[roll_indx]
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
    
    def vec_dot(a,b):
        if len(b.shape) ==3:
            c = np.zeros_like(b)
            for i in range(a.shape[0]):
                for m in range(b.shape[1]):
                    for n in range(b.shape[2]):
                        c[i,m,n] = np.sum(a[i,:]*b[:,m,n])
        elif len(b.shape) == 1:
            c = np.dot(a,b)
        return c
    
    def tp_quat2delta2(dcm_b_f,xyz,dcm_b_f_rot):
        xyz = np.array(xyz)
        vec_rot = vec_dot(dcm_b_f_rot, xyz)
        vec_f = vec_rot#vec_dot(dcm_b_f, vec_rot)###原因在于我画图这里一开始就用的统一xyz坐标，而星上坐标定义和alpha、beta不是同一坐标
        alpha = np.arctan(vec_f[0,:,:]/vec_f[2,:,:])*180/np.pi
        beta = np.arctan(vec_f[1,:,:]/vec_f[2,:,:])*180/np.pi
        return alpha,beta,vec_f[0,:,:],vec_f[1,:,:],vec_f[2,:,:]
    
    def cal_dis_meters(t1,p1,t2,p2,radius=1):
        dx = np.cos(np.pi/2.0 - t1)*np.cos(p1)-np.cos(np.pi/2.0 - t2)*np.cos(p2)
        dy = np.cos(np.pi/2.0 - t1)*np.sin(p1)-np.cos(np.pi/2.0-t2)*np.sin(p2)
        dz = np.sin(np.pi/2.0 -t1) -  np.sin(np.pi/2.0 -t2)
        del_ang = np.arcsin(np.sqrt(dx**2+dy**2+dz**2)/2)*2
        d = del_ang*radius
        return d
    
    def xyz2ab(x,y,z):
        alpha = np.arctan(x/z)
        beta = np.arctan(y/z)
        return alpha,beta
    
    def ab2xyz(alpha,beta):
        ta = np.tan(alpha)
        tb = np.tan(beta)
        z = np.sqrt(1/(ta**2+tb**2+1))
        x = z*ta
        y = z*tb
        return x, y, z
        
    def correct_psfmodel(delta_alfa0, delta_beta0, box_dex=0):
        PI = 3.14159265358979323846
        #flag = parameter[0]
        roll = 0#parameter[1] / 180.0 * PI
        roll_indx = box_dex#int(-parameter[1] / 60 + 1)
        #r_degree = parameter[1]
        psai = Tpsai[roll_indx]
        theta = Ttheta[roll_indx]
        phi = Tphi[roll_indx]
        alfa_bound = abound[roll_indx]
        beta_bound = bbound[roll_indx]
        pa, pb, pc, pd = Tpa[roll_indx], Tpb[roll_indx], Tpc[roll_indx], Tpd[roll_indx]
    
        rot_quat = quat.Quat((psai, theta, phi))
        x, y, z = ab2xyz(np.deg2rad(delta_alfa0), np.deg2rad(delta_beta0))#testroll.sph2cart(ra * pi / 180, dec * pi / 180, 1)
        xyz = [x, y, z]
        dcm_b_f_rot = rot_quat.transform.T
        delta_alfa0, delta_beta0, xr, yr, zr = tp_quat2delta2(dcm_b_f, xyz, dcm_b_f_rot)
        delta_alfa = np.abs(delta_alfa0)#np.abs(delta_alfa0 * np.cos(roll) - delta_beta0 * np.sin(roll))
        delta_beta = np.abs(delta_beta0)#np.abs(delta_alfa0 * np.sin(roll) + delta_beta0 * np.cos(roll))
        condi = np.logical_and(delta_alfa < alfa_bound, delta_beta < beta_bound)
        factor=np.where(condi, (1.0 - np.tan(delta_alfa / 180.0 * PI) / np.tan(alfa_bound / 180.0 * PI)) * \
                     (1.0 - np.tan(delta_beta / 180.0 * PI) / np.tan(beta_bound / 180.0 * PI)),np.zeros_like(delta_alfa))
        factor = factor / (np.sqrt(
            np.tan(delta_alfa / 180.0 * PI) * np.tan(delta_alfa / 180.0 * PI) + np.tan(delta_beta / 180.0 * PI) * np.tan(
                delta_beta / 180 * PI) + 1))
        rt = pa * delta_alfa ** 2 * delta_beta ** 2 + pb * delta_alfa ** 2 + pc * delta_beta ** 2 + pd
        factor = factor*rt
        return factor

    def plot_psf(func, funcname):
        fig = plt.figure(figsize=plt.figaspect(0.5))
        ax = fig.add_subplot(1, 1, 1)  # , projection='3d')
        alpha = np.arange(-Talfa_bound, Talfa_bound, 0.1)###0.001
        beta = np.arange(-Tbeta_bound, Tbeta_bound, 0.1)
        alpha, beta = np.meshgrid(alpha, beta)
        Z0 = func(alpha, beta)
        # Z1 = correct_psfmodel(alpha,beta)
        cmap = cm.coolwarm
        surf = ax.contourf(alpha, beta, Z0,cmap=cmap)  # , rstride=1, cstride=1, cmap=cm.coolwarm,
        # linewidth=0, antialiased=False)
        plt.colorbar(surf, orientation='horizontal')
        fig.savefig('{:s}_{:s}_psf.png'.format(Astr,funcname))
        return Z0

    def plot_delta(Astr):
        fig = plt.figure(figsize=plt.figaspect(0.5))
        ax = fig.add_subplot(1, 1, 1)  # , projection='3d')
        alpha = np.arange(-Talfa_bound, Talfa_bound, 0.01)
        beta = np.arange(-Tbeta_bound, Tbeta_bound, 0.01)
        alpha, beta = np.meshgrid(alpha, beta)
        Z0 = stand_psfmodel(alpha, beta)
        Z1 = correct_psfmodel(alpha, beta)
        cmap = cm.coolwarm
        ###Zratio =
        Z0 = np.ma.masked_where(Z0 == 0, Z0)
        Z1 = np.ma.masked_where(Z0 == 0, Z1)
        Zratio = Z1-Z0
        surf = ax.contourf(alpha, beta, Zratio,cmap=cmap)  # , rstride=1, cstride=1, cmap=cm.coolwarm,
        # linewidth=0, antialiased=False)
        plt.colorbar(surf, orientation='horizontal')
        fig.savefig('{:s}_psf_delta.png'.format(Astr))

    plot_psf(stand_psfmodel, 'stand')
    correct_psf = plot_psf(correct_psfmodel, 'correct')
    plot_delta(Astr)
    return correct_psf

Astr = 'A02'
correct_psf = psf(Astr)
np.save('LE_{:s}_psf'.format(Astr),correct_psf)
'''
for Astr in ['A01','A02','A03']:
    correct_psf = psf(Astr)
    np.save('{:s}_psf'.format(Astr),correct_psf)
'''

def plot_AA_delta(Z0,Z1,Astr0,Astr1):
    fig = plt.figure(figsize=plt.figaspect(0.5))
    ax = fig.add_subplot(1, 1, 1)  # , projection='3d')
    alpha = np.arange(-5.7, 5.7, 0.001)
    beta = np.arange(-1.1, 1.1, 0.001)
    alpha, beta = np.meshgrid(alpha, beta)
    cmap = cm.coolwarm
    surf = ax.contourf(alpha, beta, Z0 - Z1,cmap=cmap)  # , rstride=1, cstride=1, cmap=cm.coolwarm,
    # linewidth=0, antialiased=False)
    plt.colorbar(surf, orientation='horizontal')
    fig.savefig('delta_{:s}_{:s}.png'.format(Astr0,Astr1))
'''###
A01_crt_psf = np.load('A01_psf.npy')
A02_crt_psf = np.load('A02_psf.npy')
A03_crt_psf = np.load('A03_psf.npy')
plot_AA_delta(A01_crt_psf,A02_crt_psf,'A01','A02')
plot_AA_delta(A01_crt_psf,A03_crt_psf,'A01','A03')
plot_AA_delta(A02_crt_psf,A03_crt_psf,'A02','A03')
'''###
