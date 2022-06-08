import matplotlib.pyplot as plt

from matplotlib.transforms import (
    Bbox, TransformedBbox, blended_transform_factory)
from mpl_toolkits.axes_grid1.inset_locator import (
    BboxPatch, BboxConnector, BboxConnectorPatch)
import matplotlib.pyplot as plt
from matplotlib.pyplot import MultipleLocator,LogLocator
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


if __name__ == '__main__':

    # 原数据，共12个标记(点)
    ###Crab_ra = 83.633, Crab_dec = 22.014
    # ra = 83.633
    # dec = 22.014
    # x =  ra + np.random.uniform(low=-1.0, high=1.0, size=12)
    # y =  dec + np.random.uniform(low=-1.0, high=1.0, size=12)
    data_array = np.transpose(np.loadtxt('result_nangyi_b_5sig_for_Paper.txt'))
    exposure = data_array[0]
    static_err = data_array[1]
    sys_err = data_array[3]
    sys_bias = data_array[2]
    # dex = (x_err < 0.03) & (y_err < 0.03)
    # x = x[dex]
    # x_err = x_err[dex]
    # y = y[dex]
    # y_err = y_err[dex]

    exposure = exposure[0:11]
    static_err = static_err[0:11]
    sys_err = sys_err[0:11]
    sys_bias = sys_bias[0:11]


    # 小图扩展系数
    x_ratio = 0.8  # x轴多显示间隔
    y_ratio = 0.8  # y轴多显示间隔

    # ====== 原图 ======
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(1,1,1)
    ax.grid()


    # ax.xaxis.set_major_locator(MultipleLocator(2.5))
    # ax.yaxis.set_major_locator(MultipleLocator(2))
    ax.errorbar(  #
        exposure - 0.075*exposure, sys_bias, yerr=static_err,
        #color='blue',  #
        linestyle='', linewidth=2,  #
        marker='o', markersize=2, markeredgewidth=2,capsize=4,capthick=1,label='Statiscal error')
    ax.errorbar(  #
        exposure + 0.075*exposure, sys_bias, yerr=sys_err, #
        color='black',  #
        linestyle='', linewidth=2,  #
        marker='o', markersize=2, markeredgewidth=2, markeredgecolor='black', markerfacecolor='black',capsize=4,capthick=1,label='Systematic error')


    ax.xaxis.set_major_locator(LogLocator(base=2))
    ax.yaxis.set_major_locator(MultipleLocator(0.1))
    ax.xaxis.set_minor_locator(LogLocator(base=2))#,subs=[0.2,0.4,0.6,0.8]))
    ax.yaxis.set_minor_locator(MultipleLocator(0.02))
    ax.legend(loc='upper right', prop={'size': 12}, framealpha=1)
    ax.tick_params(axis='both', which='major', bottom=True, top=False, left=True, right=False, labelsize=14, width=2,
                   length=12)
    ax.tick_params(axis='both', which='minor', bottom=True, top=False, left=True, right=False, labelsize=7, width=1,
                   length=5)
    ax.set_xscale("log",basex=2)#,subs=[0,1,2,3,4,5,6,7,8,9,10,11,12,13])
    font2 = {
        'weight': 'normal',
        'size': 17,
    }
    # 'family': 'Times New Roman',
    ax.set_ylabel('Delta (deg)', font2)
    ax.set_xlabel('Exposure (s)', font2)
    ax.spines['top'].set_linewidth('1.2')
    ax.spines['bottom'].set_linewidth('1.2')
    ax.spines['left'].set_linewidth('1.2')
    ax.spines['right'].set_linewidth('1.2')


    plt.savefig('pointcraberror.pdf', bbox_inches='tight', dpi=fig.dpi)
    plt.show()