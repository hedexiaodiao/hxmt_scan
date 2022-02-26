import os

import matplotlib.pyplot as plt
from matplotlib.pyplot import MultipleLocator
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.transforms import TransformedBbox
from mpl_toolkits.axes_grid1.inset_locator import BboxPatch, BboxConnector
import matplotlib.patches as mpathes
import numpy as np
#  Crab, ra = 83.633, dec = 22.014
#     H 1730-333, ra = 263.353, dec = -33.3888
#     GX 354-0, ra = 262.991, dec = -33.834


if __name__ == '__main__':
    #current_dir = r'G:\ihep5\idl_psf_location\moreLE_2_6_100w'
    #current_dir = r'G:\ihep5\idl_psf_location\LE_2_6_10w'
    #current_dir = '.'
    current_dir = r'G:\ihep5\idl_psf_location\moreLE_2_6_ME_7_20_10w'
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
    burst_num = data[:,-1]

    # print(burst_num[10])
    # ax.errorbar(x[10], y[10], yerr=y_err[10], xerr=x_err[10], color='green', linewidth=1,  #
    #             linestyle='',
    #             marker='o', markersize=3, markeredgecolor='green', markerfacecolor='green')

    #================首先筛除误差过大的=============
    dex = x_err<0.5
    x = x[dex]
    x_err = x_err[dex]
    y = y[dex]
    y_err = y_err[dex]
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
    print(xG354[0], yG354[0])
    print('G354暴发序号：',G354_burst)

    # ================再筛除误差过大的=============
    dex = (x_err < 0.22) & ((x > H1730_ra - x_H1730_devia) | (y >H1730_dec - y_H1730_devia))
    x = x[dex]
    x_err = x_err[dex]
    y = y[dex]
    y_err = y_err[dex]
    burst_num = burst_num[dex]
    #==================找到合适的H1730============
    nofar_dex = (x<(H1730_ra+0.25)) & (x>G354_ra)
    x = x[nofar_dex]
    x_err = x_err[nofar_dex]
    y = y[nofar_dex]
    y_err = y_err[nofar_dex]
    burst_num = burst_num[nofar_dex]
    print('H1730暴发序号:',burst_num)
    print(x[0],y[0])

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
    ax.set_title(r"ax1")
    ax.set_xlim(H1730_ra-2, H1730_ra+2)  # 子图窗口大小固定了，调整子坐标系的显示范围
    ax.set_ylim(H1730_dec-2, H1730_dec+2)
    ax.xaxis.set_major_locator(MultipleLocator(0.5))
    ax.yaxis.set_major_locator(MultipleLocator(0.5))
    ax.plot(  #
        H1730_ra, H1730_dec,  #
        color='red',  #
        linestyle='', linewidth=2,  #
        marker='+', markersize=18, markeredgecolor='red', markerfacecolor='C0')
    ax.plot(  #
        G354_ra, G354_dec,  #
        color='red',  #
        linestyle='', linewidth=2,  #
        marker='x', markersize=18, markeredgecolor='red', markerfacecolor='C0')

    plt.show()