import matplotlib.pyplot as plt

from matplotlib.transforms import (
    Bbox, TransformedBbox, blended_transform_factory)
from mpl_toolkits.axes_grid1.inset_locator import (
    BboxPatch, BboxConnector, BboxConnectorPatch)
import matplotlib.pyplot as plt
from matplotlib.pyplot import MultipleLocator
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.transforms import TransformedBbox
from mpl_toolkits.axes_grid1.inset_locator import BboxPatch, BboxConnector
import matplotlib.patches as mpathes
import numpy as np
import sys
import matplotlib.gridspec as gridspec
import matplotlib
matplotlib.rcParams['xtick.direction'] = 'in'
matplotlib.rcParams['ytick.direction'] = 'in'


def rotate_from_mid(xy, ox, oy, theta):
    # theta must in rad
    print(xy, ox, oy)
    xy_homogeneous = np.array([xy[0], xy[1], 1])
    mtx_trans2O = np.array([[1, 0, -ox], [0, 1, -oy], [0, 0, 1]])
    mtx_rotate = np.array([[np.cos(theta), -np.sin(theta), 0], [np.sin(theta), np.cos(theta), 0], [0, 0, 1]])
    mtx_anti_trans2O = np.array([[1, 0, ox], [0, 1, oy], [0, 0, 1]])
    xy_new_homogeneous = np.dot(mtx_anti_trans2O, np.dot(mtx_rotate, np.dot(mtx_trans2O, xy_homogeneous)))
    xy_new = xy_new_homogeneous[0:2]
    return xy_new

def get_abbound(instr):
    if instr == "HE":
        Talfa_bound = 5.7  # means long_bound
        Tbeta_bound = 1.1  # means short_bound
        dcm_b_f = np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]])
    elif instr == "ME":
        Talfa_bound = 4.0  # means long_bound
        Tbeta_bound = 1.0  # means short_bound
        dcm_b_f = np.array([[0, 0, -1], [0, 1, 0], [1, 0, 0]])
    elif instr == "LE":
        Talfa_bound = 6.0  # means long_bound
        Tbeta_bound = 1.6  # means short_bound
        dcm_b_f = np.array([[0, 0, 1], [0, 1, 0], [-1, 0, 0]])
    return Talfa_bound,Tbeta_bound

def add_plot_box(ax,midx,midy,abound,bbound,roll,labelname,color,ifshow):
    xy_down_left = np.array([midx - abound, midy - bbound])

    xy_down_left = rotate_from_mid(xy_down_left, midx, midy, np.deg2rad(roll))
    xy_down_right = np.array([midx + abound, midy - bbound])
    xy_down_right = rotate_from_mid(xy_down_right, midx, midy, np.deg2rad(roll))
    xy_top_left = np.array([midx - abound, midy + bbound])
    xy_top_left = rotate_from_mid(xy_top_left, midx, midy, np.deg2rad(roll))
    xy_top_right = np.array([midx + abound, midy + bbound])
    xy_top_right = rotate_from_mid(xy_top_right, midx, midy, np.deg2rad(roll))

    x_dl_dr = np.array([xy_down_left[0], xy_down_right[0]])
    y_dl_dr = np.array([xy_down_left[1], xy_down_right[1]])
    x_dr_tr = np.array([xy_down_right[0], xy_top_right[0]])
    y_dr_tr = np.array([xy_down_right[1], xy_top_right[1]])
    x_tr_tl = np.array([xy_top_right[0], xy_top_left[0]])
    y_tr_tl = np.array([xy_top_right[1], xy_top_left[1]])
    x_tl_dl = np.array([xy_top_left[0], xy_down_left[0]])
    y_tl_dl = np.array([xy_top_left[1], xy_down_left[1]])

    if ifshow:
        rect = mpathes.Rectangle(xy_down_left, width=abound * 2, height=bbound * 2, fill=False, angle=roll, color=color, lw = 2,
                                 label=labelname)
    else:
        rect = mpathes.Rectangle(xy_down_left, width=abound * 2, height=bbound * 2, fill=False, angle=roll, color=color, lw = 2)
    ax.add_patch(rect)

if __name__ == '__main__':

    # 原数据，共12个标记(点)
    ###Crab_ra = 83.633, Crab_dec = 22.014
    ra = 83.633
    dec = 22.014
    # x =  ra + np.random.uniform(low=-1.0, high=1.0, size=12)
    # y =  dec + np.random.uniform(low=-1.0, high=1.0, size=12)
    data_array = np.transpose(np.loadtxt('Nangyi_result_le_2-6_me_7-40_he_25-100_Time_Bin-1024s_for_Paper.txt'))
    x = data_array[0]
    x_err = data_array[1]
    y = data_array[2]
    y_err = data_array[3]
    # dex = (x_err < 0.03) & (y_err < 0.03)
    # x = x[dex]
    # x_err = x_err[dex]
    # y = y[dex]
    # y_err = y_err[dex]


    # 要放大区间，索引为8-10的点
    point_first = np.argmin(x)
    point_last = np.argmin(x)

    # 小图扩展系数
    x_ratio = 0.8  # x轴多显示间隔
    y_ratio = 0.8  # y轴多显示间隔

    # --- 计算小图的坐标刻度区间 --- 【重点1】
    # 小图X轴的显示范围
    # x_diff_left = (x[point_last] - x[point_first]) * x_ratio  # 右边的比左边的大，求差值
    # x_diff_right = (x[point_last] - x[point_first]) * x_ratio  #
    # x_lim_left = x[point_first] - x_diff_left  # 左边范围往左挪半个刻度
    # x_lim_right = x[point_last] + x_diff_right  # 右边往右挪
    x_max = np.max(x)
    x_min = np.min(x)
    # y_max = max(y[point_first:(point_last+1)])
    # y_min = min(y[point_first:(point_last+1)])
    x_diff_left = (x_max - x_min) * x_ratio  # 取两个点的差
    x_diff_right = (x_max - x_min) * x_ratio  #
    x_lim_left = x_min - x_diff_left  # 底部往下挪
    x_lim_right = x_max + x_diff_right  # 顶部往上挪


    # 小图Y轴的显示范围
    y_max = np.max(y)
    y_min = np.min(y)
    # y_max = max(y[point_first:(point_last+1)])
    # y_min = min(y[point_first:(point_last+1)])
    y_diff_left = (y_max - y_min) * y_ratio  # 取两个点的差
    y_diff_right = (y_max - y_min) * y_ratio  #
    y_lim_left = y_min - y_diff_left  # 底部往下挪
    y_lim_right = y_max + y_diff_right  # 顶部往上挪

    # x_lim_left = 83.412
    # x_lim_right = 83.837
    # y_lim_left = 21.813
    # y_lim_right = 22.232
    x_lim_left = ra - 0.11
    x_lim_right = ra + 0.11
    y_lim_left = dec - 0.11
    y_lim_right = dec + 0.11

    # ====== 原图 ======
    fig = plt.figure(figsize=(12, 6))
    gs = gridspec.GridSpec(10,20)#fig.add_gridspec(10,24)
    ax = fig.add_subplot(gs[:,0:10])
    axins = fig.add_subplot(gs[1:9, 12:20])

    midx = ra
    midy = dec
    instr_list = ["HE","ME","LE"]
    roll_list = [60,0,-60]
    for instr in instr_list:
        abound, bbound = get_abbound(instr)
        if instr != "HE":
            if instr == 'ME':
                color = 'g'
                labelname = instr #+ r' 1$^\circ$.0 $\times$ 4$^\circ$.0'
            else:
                color = 'b'
                labelname = instr #+ r' 1$^\circ$.6 $\times$ 6$^\circ$.0'
            for roll in roll_list:
                ifshow = roll == 0
                roll = roll + 90
                add_plot_box(ax, midx, midy, abound, bbound, roll, labelname, color, ifshow)
        else:
            color = 'r'
            labelname = instr #+ r' 1$^\circ$.1 $\times$ 5$^\circ$.7'
            for roll in roll_list:
                ifshow = roll == 0
                # roll = roll + 90
                add_plot_box(ax, midx, midy, abound, bbound, roll, labelname, color, ifshow)


    # rect = mpathes.Rectangle(xy_down_left, width=abound*2, height=bbound*2, fill=False,angle=roll, color='r')
    # ax.add_patch(rect)

    #fig, ax = plt.subplots(1, 1, 1, figsize=(6, 4))
    # ax.set_title(r"ax1")
    ax.set_xlim(midx-8.4, midx+8.4)  # 子图窗口大小固定了，调整子坐标系的显示范围
    ax.set_ylim(midy-8.4, midy+8.4)
    ax.xaxis.set_major_locator(MultipleLocator(2.5))
    ax.yaxis.set_major_locator(MultipleLocator(2))
    ax.errorbar(  #
        x, y, xerr=x_err, yerr=y_err, #
        color='black',  #
        linestyle='', linewidth=1,  #
        marker='o', markersize=2, markeredgecolor='black', markerfacecolor='black')
    ax.plot(  #
        ra, dec,
        color='red',  #
        linestyle='', linewidth=3,  #
        marker='+', markersize=8, markeredgewidth=2, markeredgecolor='red', markerfacecolor='red',label='Crab')

    # ====== 子图 ======
    # 插入子图。绘制局部放大图的坐标系
    # axins = inset_axes(  #
    #     ax,  #
    #     width="40%", height="30%",  # 宽高，str百分比或float小数比例
    #     loc=1,  # 子图锚点位置：
    #     # 'upper right'  : 1, 右上角
    #     # 'upper left'   : 2, 左上角
    #     # 'lower left'   : 3, 左下角
    #     # 'lower right'  : 4, 右下角
    #     # 'center left'  : 6, 左中
    #     # 'center right' : 7, 右中
    #     # 'lower center' : 8, 底中
    #     # 'upper center' : 9, 顶中
    #     # 'center'       : 10 图中心
    #     bbox_to_anchor=(0, -0.1, 1, 1),  # 根据锚点设置子图的位置
    #     # loc设置的锚点系，决定了这里的便宜方向 。
    #     # left, 锚点向【右】偏离原位置的距离比例，窗口宽的百分比 0-1
    #     # bottom, 锚点向【上】偏离原位置的距离比例，窗口高的百分比 0-1
    #     # width, [保持1]，宽使用参数width的设置
    #     # height, [保持1]，高使用参数height的设置
    #     bbox_transform=ax.transAxes)


    axins.errorbar(  # 子图绘制原始数据
        x, y, xerr=x_err, yerr=y_err,zorder = 1, #
        color='black', linestyle='', linewidth=1,  #
        marker='o', markersize=2, markeredgecolor='black', markerfacecolor='black',label='Fitting Positions')
    # axins.set_title("ax2")
    axins.plot(  #
        ra, dec, zorder = 2,
        color='red',  #
        linestyle='', linewidth=3,  #
        marker='+', markersize=16, markeredgewidth=2, markeredgecolor='red', markerfacecolor='red',label='Crab')

    # --- 设置小图的坐标刻度区间 --- 【重点1】
    print(x_lim_left,x_lim_right)
    axins.xaxis.set_major_locator(MultipleLocator(0.04))
    axins.yaxis.set_major_locator(MultipleLocator(0.04))
    axins.set_xlim(x_lim_left, x_lim_right)  # 子图窗口大小固定了，调整子坐标系的显示范围
    axins.set_ylim(y_lim_left, y_lim_right)
    ###axins.xaxis.set_major_locator(MultipleLocator(0.25))
    ###axins.yaxis.set_major_locator(MultipleLocator(0.01))
    ax.legend(loc='upper left', prop={'size': 12}, framealpha=1)
    axins.legend(loc='upper left', prop={'size': 12}, framealpha=1)
    ax.tick_params(axis='both', which='major', bottom=True, top=False, left=True, right=False, labelsize=14, width=2,
                   length=10)
    axins.tick_params(axis='both', which='major', bottom=True, top=False, left=True, right=False, labelsize=14, width=2,
                   length=10)
    ax.spines['top'].set_linewidth('1.2')
    ax.spines['bottom'].set_linewidth('1.2')
    ax.spines['left'].set_linewidth('1.2')
    ax.spines['right'].set_linewidth('1.2')
    axins.spines['top'].set_linewidth('1.2')
    axins.spines['bottom'].set_linewidth('1.2')
    axins.spines['left'].set_linewidth('1.2')
    axins.spines['right'].set_linewidth('1.2')


    # 设置横纵坐标的名称以及对应字体格式
    font2 = {
             'weight': 'normal',
             'size': 17,
             }
    #'family': 'Times New Roman',
    ax.set_xlabel('RA (deg)',font2)
    ax.set_ylabel('DEC (deg)', font2)
    axins.set_xlabel('RA (deg)', font2)
    axins.set_ylabel('DEC (deg)', font2)

    # ====== 原图的要放大区域 ======
    rect = TransformedBbox(axins.viewLim, ax.transData)
    box = BboxPatch(rect, color="r", alpha=0.0, ec="b", lw=5)

    ax.add_patch(box)

    # ====== 2个图的连接线 ====== 【重点2】-BboxConnector
    line_1 = BboxConnector(  # 连线
        axins.bbox,  # bbox1
        rect,  # bbox2
        loc1=2,  # 从bbox1的哪个点连线。右上:1, 左上:2, 左下:3, 右下:4
        loc2=1,  # 从bbox2的哪个点连线。
        color='g',
        lw=3
    )
    line_1.set_clip_on(False)  # 关闭图形重叠时的布尔运算。否则线条会不显示
    line_2 = BboxConnector(  # 连线
        axins.bbox,  # bbox1
        rect,  # bbox2
        loc1=3,  # 从bbox1的哪个点连线。右上:1, 左上:2, 左下:3, 右下:4
        loc2=4,  # 从bbox2的哪个点连线。
        color='g',
        lw=3
    )
    line_2.set_clip_on(False)

    axins.add_patch(line_1)
    axins.add_patch(line_2)

    plt.savefig('Crab1024s.pdf', bbox_inches='tight', dpi=fig.dpi)
    plt.show()