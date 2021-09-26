#!/home/hxmt/nangyi/anaconda2/bin/bin/python
from xspec import *
import math
import numpy as np
import time
from astropy.io import fits as pf
from astropy import units as u
from astropy.coordinates import SkyCoord
import sys,os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
matplotlib.rcParams['xtick.direction'] = 'in'
matplotlib.rcParams['ytick.direction'] = 'in'
#from sympy import *
###########Class and Function defined by HDPC################################


def dataplot1(d1,d2,d3,m1,m2,m3,infilePri):
    fig=plt.figure(figsize=(20,10))
    gs = gridspec.GridSpec(3, 1,bottom=0.05,top=0.96,left=0.05,right=0.95,hspace=0.05)
    gs0 = gridspec.GridSpecFromSubplotSpec(2, 1,subplot_spec=gs[0], hspace=0.)
    ax00 = plt.subplot(gs0[0])
    ax01 = plt.subplot(gs0[1],sharex=ax00)
    gs1 = gridspec.GridSpecFromSubplotSpec(2, 1,subplot_spec=gs[1], hspace=0.)
    ax10 = plt.subplot(gs1[0],sharex=ax00)
    ax11 = plt.subplot(gs1[1],sharex=ax00)
    gs2 = gridspec.GridSpecFromSubplotSpec(2, 1,subplot_spec=gs[2], hspace=0.)
    ax20 = plt.subplot(gs2[0],sharex=ax00)
    ax21 = plt.subplot(gs2[1],sharex=ax00)
    ax00.plot(d1.noticed, np.array(d1.values),'b-',lw=1,label='data')
    ax10.plot(d2.noticed, np.array(d2.values),'b-',lw=1)
    ax20.plot(d3.noticed, np.array(d3.values),'b-',lw=1)
    ax00.plot(d1.noticed, np.array(m1.folded(1)),'r-',lw=2.5,label='PSF')
    ax10.plot(d2.noticed, np.array(m2.folded(2)),'r-',lw=2.5)
    ax20.plot(d3.noticed, np.array(m3.folded(3)),'r-',lw=2.5)
    ax01.plot(d1.noticed, np.array(d1.values)- np.array(m1.folded(1)),'k-')
    ax11.plot(d2.noticed, np.array(d2.values)- np.array(m2.folded(2)),'k-')
    ax21.plot(d3.noticed, np.array(d3.values)- np.array(m3.folded(3)),'k-')
    lgd=ax00.legend(loc='upper right',fontsize=13)
    ax00.set_title(infilePri+" Count Rates Curve",fontsize=15)
    ax00.set_ylabel('Grp0 Rates(cts/s)',fontsize=10)
    ax10.set_ylabel('Grp1 Rates(cts/s)',fontsize=10)
    ax20.set_ylabel('Grp5 Rates(cts/s)',fontsize=10)
    ax01.set_ylabel('Residuals',fontsize=10)
    ax11.set_ylabel('Residuals',fontsize=10)
    ax21.set_ylabel('Residuals',fontsize=10)
    ax21.set_xlabel('Times(s)',fontsize=18)
    ax00.tick_params(labelbottom=False)
    ax01.tick_params(labelbottom=False)
    ax11.tick_params(labelbottom=False)
    ax10.tick_params(labelbottom=False)
    ax20.tick_params(labelbottom=False)
    ax01.get_yticklabels()[-1].set_color('None')
    ax11.get_yticklabels()[-1].set_color('None')
    ax21.get_yticklabels()[-1].set_color('None')
    plt.xlim(0,len(d1.noticed))
    ax00.grid()
    ax10.grid()
    ax20.grid()
    ax01.grid()
    ax11.grid()
    ax21.grid()
    axe = [ax00,ax10,ax20]
    max1 = max(np.array(d1.values))
    max2 = max(np.array(d2.values))
    max3 = max(np.array(d3.values))
    min1 = min(np.array(d1.values))
    min2 = min(np.array(d2.values))
    min3 = min(np.array(d3.values))
    max1 = max1 + 18
    max2 = max2 + 18
    max3 = max3 + 18
    maxt = [max1,max2,max3]
    axe[0].set_ylim(min1-5,max1)
    axe[1].set_ylim(min2-5,max2)
    axe[2].set_ylim(min3-5,max3)
    plt.savefig('%s_KnownLc.eps'%(infilePri),dpi = fig.dpi)
    os.system('convert ' + '%s_KnownLc.eps'%(infilePri) + ' ' + '%s_KnownLc.png'%(infilePri))
    plt.close('all')

def dataplot2(d1,d2,d3,m1,m2,m3,contri1,contri2,contri3,infilePri,i,labelnm):
    fig = plt.figure(figsize=(20, 10))
    gs = gridspec.GridSpec(3, 1, bottom=0.05, top=0.96, left=0.05, right=0.95, hspace=0.05)
    gs0 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs[0], hspace=0.)
    ax00 = plt.subplot(gs0[0])
    ax01 = plt.subplot(gs0[1], sharex=ax00)
    gs1 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs[1], hspace=0.)
    ax10 = plt.subplot(gs1[0], sharex=ax00)
    ax11 = plt.subplot(gs1[1], sharex=ax00)
    gs2 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs[2], hspace=0.)
    ax20 = plt.subplot(gs2[0], sharex=ax00)
    ax21 = plt.subplot(gs2[1], sharex=ax00)
    # fig,axe = plt.subplots(3,sharex=True, sharey=False,figsize=(20,16))

    ax00.plot(d1.noticed, np.array(d1.values), 'b-', lw=1, label='data')
    ax10.plot(d2.noticed, np.array(d2.values), 'b-', lw=1)
    ax20.plot(d3.noticed, np.array(d3.values), 'b-', lw=1)
    ax00.plot(d1.noticed, np.array(m1.folded(1)), 'k-.', lw=2, label='All Src')
    ax10.plot(d2.noticed, np.array(m2.folded(2)), 'k-.', lw=2)
    ax20.plot(d3.noticed, np.array(m3.folded(3)), 'k-.', lw=2)
    ax01.plot(d1.noticed, np.array(d1.values) - contri1 - m1(2).values[0], 'k-')
    ax11.plot(d2.noticed, np.array(d2.values) - contri2 - m2(2).values[0], 'k-')
    ax21.plot(d3.noticed, np.array(d3.values) - contri3 - m3(2).values[0], 'k-')
    ax00.plot(d1.noticed, contri1 + m1(2).values[0], 'r-', lw=2.5, label='%s' % labelnm)
    ax10.plot(d2.noticed, contri2 + m2(2).values[0], 'r-', lw=2.5)
    ax20.plot(d3.noticed, contri3 + m3(2).values[0], 'r-', lw=2.5)
    lgd = ax00.legend(loc='upper right', fontsize=13)
    ax00.set_title(infilePri + "_%s Count Rates Curve" % (int((i - 2) / 5)), fontsize=15)
    ax00.set_ylabel('Grp0 Rates(cts/s)', fontsize=10)
    ax10.set_ylabel('Grp1 Rates(cts/s)', fontsize=10)
    ax20.set_ylabel('Grp5 Rates(cts/s)', fontsize=10)
    ax01.set_ylabel('Residuals', fontsize=10)
    ax11.set_ylabel('Residuals', fontsize=10)
    ax21.set_ylabel('Residuals', fontsize=10)
    ax21.set_xlabel('Times(s)', fontsize=18)
    ax00.tick_params(labelbottom=False)
    ax01.tick_params(labelbottom=False)
    ax11.tick_params(labelbottom=False)
    ax10.tick_params(labelbottom=False)
    ax20.tick_params(labelbottom=False)
    ax01.get_yticklabels()[-1].set_color('None')
    ax11.get_yticklabels()[-1].set_color('None')
    ax21.get_yticklabels()[-1].set_color('None')
    plt.xlim(0, len(d1.noticed))
    ax00.grid()
    ax10.grid()
    ax20.grid()
    ax01.grid()
    ax11.grid()
    ax21.grid()
    frame = lgd.get_frame()
    frame.set_alpha(0.5)
    # fig.subplots_adjust(hspace= 0.05)

    plt.savefig('%s_%s_src.eps' % (infilePri, int((i - 2) / 5)), dpi=fig.dpi)
    os.system('convert ' + '%s_%s_src.eps' % (infilePri, int((i - 2) / 5)) + ' ' + '%s_%s_src.png' % (infilePri, int((i - 2) / 5)))
    plt.close('all')

def dataplot3(d1,d2,d3,m1,m2,m3,infilePri,dtemp_list,srcpot,srclist):
    fig = plt.figure(figsize=(20, 10))
    gs = gridspec.GridSpec(3, 1, bottom=0.05, top=0.96, left=0.05, right=0.95, hspace=0.05)
    gs0 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs[0], hspace=0.)
    ax00 = plt.subplot(gs0[0])
    ax01 = plt.subplot(gs0[1], sharex=ax00)
    gs1 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs[1], hspace=0.)
    ax10 = plt.subplot(gs1[0], sharex=ax00)
    ax11 = plt.subplot(gs1[1], sharex=ax00)
    gs2 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs[2], hspace=0.)
    ax20 = plt.subplot(gs2[0], sharex=ax00)
    ax21 = plt.subplot(gs2[1], sharex=ax00)
    ax00.plot(d1.noticed, np.array(d1.values), 'b-', lw=1, label='data')
    ax10.plot(d2.noticed, np.array(d2.values), 'b-', lw=1)
    ax20.plot(d3.noticed, np.array(d3.values), 'b-', lw=1)
    ax00.plot(d1.noticed, np.array(m1.folded(1)), 'r-', lw=2.5, label='PSF')
    ax10.plot(d2.noticed, np.array(m2.folded(2)), 'r-', lw=2.5)
    ax20.plot(d3.noticed, np.array(m3.folded(3)), 'r-', lw=2.5)
    ax01.plot(d1.noticed, np.array(d1.values) - np.array(m1.folded(1)), 'k-')
    ax11.plot(d2.noticed, np.array(d2.values) - np.array(m2.folded(2)), 'k-')
    ax21.plot(d3.noticed, np.array(d3.values) - np.array(m3.folded(3)), 'k-')
    lgd = ax00.legend(loc='upper right', fontsize=13)
    ax00.set_title(infilePri + " Count Rates Curve", fontsize=15)
    ax00.set_ylabel('Grp0 Rates(cts/s)', fontsize=10)
    ax10.set_ylabel('Grp1 Rates(cts/s)', fontsize=10)
    ax20.set_ylabel('Grp5 Rates(cts/s)', fontsize=10)
    ax01.set_ylabel('Residuals', fontsize=10)
    ax11.set_ylabel('Residuals', fontsize=10)
    ax21.set_ylabel('Residuals', fontsize=10)
    ax21.set_xlabel('Times(s)', fontsize=18)
    ax00.tick_params(labelbottom=False)
    ax01.tick_params(labelbottom=False)
    ax11.tick_params(labelbottom=False)
    ax10.tick_params(labelbottom=False)
    ax20.tick_params(labelbottom=False)
    ax01.get_yticklabels()[-1].set_color('None')
    ax11.get_yticklabels()[-1].set_color('None')
    ax21.get_yticklabels()[-1].set_color('None')
    plt.xlim(0, len(d1.noticed))
    ax00.grid()
    ax10.grid()
    ax20.grid()
    ax01.grid()
    ax11.grid()
    ax21.grid()
    axe = [ax00, ax10, ax20]
    max1 = max(np.array(d1.values))
    max2 = max(np.array(d2.values))
    max3 = max(np.array(d3.values))
    min1 = min(np.array(d1.values))
    min2 = min(np.array(d2.values))
    min3 = min(np.array(d3.values))
    max1 = max1 + 18
    max2 = max2 + 18
    max3 = max3 + 18
    maxt = [max1, max2, max3]
    axe[0].set_ylim(min1 - 5, max1)
    axe[1].set_ylim(min2 - 5, max2)
    axe[2].set_ylim(min3 - 5, max3)

    ratpanot = 0
    anotx = [0, 9999999999]
    anotxr = []
    anote = [0]
    anoter = []
    tplen = (len(srcpot) / 3)
    for j in range(3):
        for i in range(int(tplen)):
            for h, potp in enumerate(srcpot[3 * i + j]):
                for index, l in enumerate(anotx):
                    if potp <= l:
                        anotx.insert(index, potp)
                        anote.insert(index, i + 1)
                        break
        anotxr.append(anotx[1:-1])
        anoter.append(anote[1:])
        anotx = [0, 9999999999]
        anote = [0]
    # high = [max1-17,max2-17,max3-17]
    for j in range(3):
        anotlen = len(anotxr[j])
        if anotlen < 24:
            locate = np.linspace(5, len(d3.noticed), 25)[:-1]
        else:
            locate = np.linspace(5, len(d3.noticed), anotlen + 1)[:-1]
        # if len(anotxr[j])>12:
        #    tmlt = range(len(anotxr[j]))
        #    for tm in range(len(anotxr[j])):
        #        tmlt[tm] = high[j] if ((tm+1) % 2)==0 else 0
        #        print j, "--",tm
        # else:
        #    tmlt = [high[j] for tm in range(len(anotxr[j]))]
        for i, potp in enumerate(anotxr[j]):
            dtemp = dtemp_list[j]
            hpot = (np.array(dtemp.folded(j + 1))[(np.array(d3.noticed) == (int(potp)))])[0]
            locnow = locate[np.argmin(np.abs(locate - potp))]
            locate = np.delete(locate, np.argmin(np.abs(locate - potp)))
            axe[j].annotate("(%s)"%((srclist[anoter[j][i]-1]-2)/5),xy=(potp,hpot),xytext=(locnow,maxt[j]*7./10.0),textcoords="data",fontsize=16,arrowprops=dict(arrowstyle='->'),color='k')

    frame = lgd.get_frame()
    frame.set_alpha(0.4)
    plt.savefig('%s_LcFit.eps' % (infilePri), dpi=fig.dpi)
    os.system('convert ' + '%s_LcFit.eps' % (infilePri) + ' ' + '%s_LcFit.png' % (infilePri))
    plt.close('all')

def dataplot4(infilePri,srcs_bright,ratp,dectp,goodra,gooddec,goodsnr,goodnorm,src_map_near,time1):
    cm = plt.cm.get_cmap('RdBu_r')  ###
    fig, axe = plt.subplots(1, 1, figsize=(24, 15))  ###
    axe.set_facecolor('gray')
    if len(srcs_bright[:, 1]) > 0:
        plt.scatter(srcs_bright[:, 1], srcs_bright[:, 2], c=srcs_bright[:, 0], marker='o', label="Source File bright",
                    cmap=cm, s=60)

    axe.plot(ratp, dectp, 'b>', alpha=0.6, ms=0.5)
    clb = plt.scatter(np.array(goodra)[(np.array(goodsnr) < 10000)], np.array(gooddec)[(np.array(goodsnr) < 10000)],
                      c=np.array(goodnorm)[(np.array(goodsnr) < 10000)], marker='*', label="snr>5", s=300, cmap=cm,
                      alpha=.5)
    cbar = plt.colorbar(clb, pad=0.01)
    ###'''
    time2 = time.time()
    with open("process", 'a') as f:
        f.write("Test! 2\t---> %s" % ((time2 - time1) / 3600.0) + ' hours' + str(
            time.strftime('%H:%M.%S, %Y.%m.%d', time.localtime(time.time()))) + "\n")

    ###'''
    for i, src in enumerate(srcs_bright):
        plt.text(src[1] + 0.1, src[2], "%s" % (i + 1), fontsize=15)

    plt.scatter(src_map_near[:, 1], src_map_near[:, 2], c='w', marker='.', label="Source File", s=15)  # cmap=cm
    ###'''

    ###'''
    axe.set_title(infilePri, fontsize=12)
    plt.legend()
    plt.grid()
    cbar.set_label('Norm', fontsize=15)
    axe.set_xlabel('Ra', fontsize=20)
    axe.set_ylabel('Dec', fontsize=20)
    plt.tight_layout()
    plt.savefig('%sA-S.eps' % infilePri, dpi=fig.dpi)
    os.system('convert '+'%sA-S.eps' % infilePri+' '+'%sA-S.png' % infilePri)
    # plt.show()
    time2 = time.time()
    with open("process", 'a') as f:
        f.write("Test! 4\t---> %s" % ((time2 - time1) / 3600.0) + ' hours' + str(
            time.strftime('%H:%M.%S, %Y.%m.%d', time.localtime(time.time()))) + "\n")

    plt.close('all')
