import os

import matplotlib.pyplot as plt
from matplotlib.pyplot import MultipleLocator
import matplotlib
matplotlib.rcParams['xtick.direction'] = 'in'
matplotlib.rcParams['ytick.direction'] = 'in'
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.transforms import TransformedBbox
from mpl_toolkits.axes_grid1.inset_locator import BboxPatch, BboxConnector
import matplotlib.patches as mpathes
import numpy as np
#  Crab, ra = 83.633, dec = 22.014
#     H 1730-333, ra = 263.353, dec = -33.3888
#     GX 354-0, ra = 262.991, dec = -33.834
def cal_dis_meters(t1, p1, t2, p2, radius=1):
    t1 = np.deg2rad(t1)
    p1 = np.deg2rad(p1)
    t2 = np.deg2rad(t2)
    p2 = np.deg2rad(p2)
    dx = np.cos(np.pi/2.0 - t1)*np.cos(p1)-np.cos(np.pi/2.0-t2)*np.cos(p2)#checked,OK
    dy = np.cos(np.pi/2.0 - t1)*np.sin(p1)-np.cos(np.pi/2.0-t2)*np.sin(p2)
    dz = np.sin(np.pi/2.0 - t1) - np.sin(np.pi/2.0-t2)
    del_ang = np.rad2deg(np.arcsin(np.sqrt(dx**2+dy**2+dz**2)/2)*2)
    d = del_ang*radius
    return del_ang,d

if __name__ == '__main__':
    #current_dir = r'G:\ihep5\idl_psf_location\moreLE_2_6_100w'
    #current_dir = r'G:\ihep5\idl_psf_location\LE_2_6_10w'
    #current_dir = '.'
    #current_dir = r'G:\ihep5\idl_psf_location\moreLE_2_6_ME_7_20_10w'
    #current_dir = r'G:\ihep5\idl_psf_location\moreLE_2_6_ME_7_20_10w'
    #current_dir = r'G:\ihep5\idl_psf_location\moreLE_2_6_ME_7_20_10w'
    #current_dir = r'G:\ihep5\idl_psf_location\ME_7_20_100w'
    current_dir = r'G:\ihep5\idl_psf_location\bf_v7\right_LE_2-6_ME_7-20_10w'
    #current_dir = r'G:\ihep5\idl_psf_location\bf_v7\right_LE_2-6_10w'
    #current_dir = r'G:\ihep5\idl_psf_location\right_LE_2-6_ME_7-20_10w_24b100w'
    #current_dir = r'G:\ihep5\idl_psf_location\v7_LE_2-6_ME_7-12_10w'
    #current_dir = r'G:\ihep5\idl_psf_location\v7_LE_2-6_ME_7-20_10w'
    #current_dir = r'G:\ihep5\idl_psf_location\bf_v7\right_LE_2-6_ME_7-12_10w'
    #current_dir = r'G:\ihep5\idl_psf_locnew\LE_2-6_ME_7-20_10w'
    #current_dir = r'G:\ihep5\idl_psf_locnew\LE_2-6_ME_7-12_10w'
    current_dir = r'G:\ihep5\idl_psf_locnew\LE_2-6_ME_7-40_10w'
    os.chdir(current_dir)
    #txtname = 'plotall_10w_LE_2-6_ME_7-20.txt'
    #txtname = 'plotall_LE_2-6_ME_7-20.txt'
    #txtname = 'plotall_ME_7-20.txt'
    txtname = 'plotall.txt'
    # 原数据，共12个标记(点)
    H1730_ra = 263.353
    H1730_dec = -33.3888
    G354_ra = 262.991
    G354_dec = -33.834

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(1, 1, 1)

    data =  np.loadtxt(txtname)#H1730_ra + np.random.uniform(low=-1.0, high=1.0, size=12)
    print(data.shape)
    x = data[:,0]
    x_err = data[:,1]
    y =  data[:,2]#H1730_dec + np.random.uniform(low=-1.0, high=1.0, size=12)
    y_err = data[:,3]
    deg_err = data[:,4]
    burst_num = data[:,-1]

    ax.plot(  #
        H1730_ra, H1730_dec,  #
        color='red',  #
        linestyle='', linewidth=2,  #
        marker='D', markersize=7, markeredgecolor='red', markerfacecolor='red', label='H1730-333')
    ax.plot(  #
        G354_ra, G354_dec,  #
        color='red',  #
        linestyle='', linewidth=2,  #
        marker='s', markersize=7, markeredgecolor='red', markerfacecolor='red', label='G354-0')

    # print(burst_num[10])
    # ax.errorbar(x[10], y[10], yerr=y_err[10], xerr=x_err[10], color='green', linewidth=1,  #
    #             linestyle='',
    #             marker='o', markersize=3, markeredgecolor='green', markerfacecolor='green')

    #================首先筛除误差过大的=============
    ###dex = (x_err<0.2) & (y_err<0.6) & (x>261.5) & (x<265.0) & (y>-35.0) & (y<-32.0)
    dex = (x_err < 0.23) & (y_err < 0.18) & (x > 261.5) & (x < 265.0) & (y > -35.0) & (y < -32.0)
    x = x[dex]
    x_err = x_err[dex]
    y = y[dex]
    y_err = y_err[dex]
    deg_err = deg_err[dex]
    burst_num = burst_num[dex]

    #==================找到G354==================
    x_H1730_devia = 0.1
    y_H1730_devia = 0.3
    G354_dex = (x < H1730_ra - x_H1730_devia) & (y<H1730_dec-y_H1730_devia)#(x < G354_ra)
    G354_burst = burst_num[G354_dex]
    xG354 = x[G354_dex]
    x_errG354 = x_err[G354_dex]
    yG354 = y[G354_dex]
    y_errG354 = y_err[G354_dex]
    deg_errG354 = deg_err[G354_dex]
    print(xG354[0], yG354[0])
    print(x_errG354,y_errG354)
    print('G354暴发序号：',G354_burst)
    print('G354对应ra：',xG354)
    print('G354对应ra_err：', x_errG354)
    print('G354对应dec:',yG354)
    print('G354对应dec_err:',y_errG354)
    print('G354对应序号误差：',deg_errG354)
    del_ang,d = cal_dis_meters(90 - G354_dec, G354_ra, 90 - yG354, xG354)
    print('G354对应序号偏离：', del_ang)
    # ================再筛除误差过大的=============
    dex = (x_err < 0.23) & (y_err < 0.23) & ((x > H1730_ra - x_H1730_devia) | (y >H1730_dec - y_H1730_devia))
    x = x[dex]
    x_err = x_err[dex]
    y = y[dex]
    y_err = y_err[dex]
    deg_err = deg_err[dex]
    burst_num = burst_num[dex]
    #==================找到合适的H1730============
    nofar_dex = (x<(H1730_ra+0.25)) & (x>G354_ra)
    x = x[nofar_dex]
    x_err = x_err[nofar_dex]
    y = y[nofar_dex]
    y_err = y_err[nofar_dex]
    deg_err = deg_err[nofar_dex]
    burst_num = burst_num[nofar_dex]
    print('H1730暴发序号:',burst_num)
    print('H1730对应暴发ra:',x)
    print('H1730对应暴发ra_err:',x_err)
    print('H1730对应暴发dec:',y)
    print('H1730对应暴发dec_err:', y_err)
    print('H1730对应序号误差：',deg_err)
    print(x[0],y[0])
    del_ang, d = cal_dis_meters(90 - H1730_dec, H1730_ra, 90 - y, x)
    print('H1730对应序号偏离：', del_ang)
    print('H1730对应序号平均偏离：', np.mean(del_ang))
    del_ang, d = cal_dis_meters(90 - H1730_dec, H1730_ra, 90 - G354_dec, G354_ra)
    print('H1730与G354距离：', del_ang)

    # # 要放大区间，索引为8-10的点
    # point_first = 0
    # point_last = 11
    #
    # # 小图扩展系数
    # x_ratio = 0.3  # x轴多显示间隔
    # y_ratio = 0.5  # y轴多显示间隔
    #
    # # --- 计算小图的坐标刻度区间 --- 【重点1】
    # # 小图X轴的显示范围
    # x_diff_left = (x[point_last] - x[point_first]) * x_ratio  # 右边的比左边的大，求差值
    # x_diff_right = (x[point_last] - x[point_first]) * x_ratio  #
    # x_lim_left = x[point_first] - x_diff_left  # 左边范围往左挪半个刻度
    # x_lim_right = x[point_last] + x_diff_right  # 右边往右挪
    #
    # # 小图Y轴的显示范围
    # y_max = max(y[point_first:(point_last+1)])
    # y_min = min(y[point_first:(point_last+1)])
    # y_diff_left = (y_max - y_min) * y_ratio  # 取两个点的差
    # y_diff_right = (y_max - y_min) * y_ratio  #
    # y_lim_left = y_min - y_diff_left  # 底部往下挪
    # y_lim_right = y_max + y_diff_right  # 顶部往上挪

    # ====== 原图 ======


    ax.errorbar(x, y, yerr=y_err,xerr=x_err, color='black', linewidth=1,  #
                linestyle='',
            marker='o', markersize=3, markeredgecolor='black', markerfacecolor='black')

    ax.errorbar(xG354, yG354, yerr=y_errG354, xerr=x_errG354, color='blue', linewidth=1,  #
                linestyle='',
                marker='o', markersize=3, markeredgecolor='blue', markerfacecolor='blue')

    # rect = mpathes.Rectangle(xy_down_left, width=abound*2, height=bbound*2, fill=False,angle=roll, color='r')
    # ax.add_patch(rect)

    #fig, ax = plt.subplots(1, 1, 1, figsize=(6, 4))
    #ax.set_title(r"ax1")
    ax.set_xlim(H1730_ra-2, H1730_ra+2)  # 子图窗口大小固定了，调整子坐标系的显示范围
    ax.set_ylim(H1730_dec-2, H1730_dec+2)
    ax.xaxis.set_major_locator(MultipleLocator(1.0))
    ax.yaxis.set_major_locator(MultipleLocator(0.5))
    ax.legend(loc='upper right',prop={'size':16},framealpha=1)
    ax.tick_params(axis='both', which='major', bottom=True, top=False, left=True, right=False, labelsize=14, width=1,
                   length=10)
    ax.set_xlabel('RA (deg)', fontsize=16)  # 28
    ax.set_ylabel('DEC (deg)', fontsize=16)
    plt.savefig('burst_points.pdf', bbox_inches='tight', dpi=fig.dpi)
    plt.show()