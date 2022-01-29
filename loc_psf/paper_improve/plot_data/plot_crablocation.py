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

def add_plot_box(ax,midx,midy,abound,bbound,roll):
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

    ax.plot(midx, midy, color='k', linewidth=1,  #
            marker='o', markersize=5, markeredgecolor='black', markerfacecolor='C0')
    ax.plot(x_dl_dr, y_dl_dr, color='k', linewidth=1,  #
            marker='o', markersize=5, markeredgecolor='black', markerfacecolor='C0')
    ax.plot(x_dr_tr, y_dr_tr, color='k', linewidth=1,  #
            marker='o', markersize=5, markeredgecolor='black', markerfacecolor='C0')
    ax.plot(x_tr_tl, y_tr_tl, color='k', linewidth=1,  #
            marker='o', markersize=5, markeredgecolor='black', markerfacecolor='C0')
    ax.plot(x_tl_dl, y_tl_dl, color='k', linewidth=1,  #
            marker='o', markersize=5, markeredgecolor='black', markerfacecolor='C0')

if __name__ == '__main__':

    # 原数据，共12个标记(点)
    ra = 83.633
    dec = 22.014
    x =  ra + np.random.uniform(low=-1.0, high=1.0, size=12)
    y =  dec + np.random.uniform(low=-1.0, high=1.0, size=12)

    # 要放大区间，索引为8-10的点
    point_first = 0
    point_last = 11

    # 小图扩展系数
    x_ratio = 0.3  # x轴多显示间隔
    y_ratio = 0.5  # y轴多显示间隔

    # --- 计算小图的坐标刻度区间 --- 【重点1】
    # 小图X轴的显示范围
    x_diff_left = (x[point_last] - x[point_first]) * x_ratio  # 右边的比左边的大，求差值
    x_diff_right = (x[point_last] - x[point_first]) * x_ratio  #
    x_lim_left = x[point_first] - x_diff_left  # 左边范围往左挪半个刻度
    x_lim_right = x[point_last] + x_diff_right  # 右边往右挪

    # 小图Y轴的显示范围
    y_max = max(y[point_first:(point_last+1)])
    y_min = min(y[point_first:(point_last+1)])
    y_diff_left = (y_max - y_min) * y_ratio  # 取两个点的差
    y_diff_right = (y_max - y_min) * y_ratio  #
    y_lim_left = y_min - y_diff_left  # 底部往下挪
    y_lim_right = y_max + y_diff_right  # 顶部往上挪

    # ====== 原图 ======
    fig = plt.figure(figsize=(16, 8))
    gs = fig.add_gridspec(10,24)
    ax = fig.add_subplot(gs[:,0:12])
    axins = fig.add_subplot(gs[1:9, 14:23])

    midx = ra
    midy = dec
    instr_list = ["HE","ME","LE"]
    roll_list = [-60,0,60]
    for instr in instr_list:
        abound, bbound = get_abbound(instr)
        if instr != "HE":
            for roll in roll_list:
                roll = roll + 90
                add_plot_box(ax, midx, midy, abound, bbound, roll)
        else:
            for roll in roll_list:
                #roll = roll + 90
                add_plot_box(ax, midx, midy, abound, bbound, roll)


    # rect = mpathes.Rectangle(xy_down_left, width=abound*2, height=bbound*2, fill=False,angle=roll, color='r')
    # ax.add_patch(rect)

    #fig, ax = plt.subplots(1, 1, 1, figsize=(6, 4))
    ax.set_title(r"ax1")
    ax.set_xlim(midx-8, midx+8)  # 子图窗口大小固定了，调整子坐标系的显示范围
    ax.set_ylim(midy-8, midy+8)
    ax.xaxis.set_major_locator(MultipleLocator(2.5))
    ax.yaxis.set_major_locator(MultipleLocator(2))
    ax.plot(  #
        x, y,  #
        color='k',  #
        linestyle=':', linewidth=1,  #
        marker='o', markersize=5, markeredgecolor='black', markerfacecolor='C0')

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

    axins.plot(  # 子图绘制原始数据
        x, y,  #
        color='k', linestyle=':', linewidth=1,  #
        marker='o', markersize=5, markeredgecolor='black', markerfacecolor='C0')
    axins.set_title("ax2")

    # --- 设置小图的坐标刻度区间 --- 【重点1】
    axins.set_xlim(x_lim_left, x_lim_right)  # 子图窗口大小固定了，调整子坐标系的显示范围
    axins.set_ylim(y_lim_left, y_lim_right)
    ###axins.xaxis.set_major_locator(MultipleLocator(0.25))
    ###axins.yaxis.set_major_locator(MultipleLocator(0.01))

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

    plt.show()