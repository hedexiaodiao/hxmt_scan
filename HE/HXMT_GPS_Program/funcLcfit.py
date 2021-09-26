#!/home/hxmt/nangyi/anaconda2/bin/bin/python
from xspec import *
import math
import numpy as np
import time
from astropy.io import fits as pf
from astropy import units as u
from astropy.coordinates import SkyCoord
import sys, os
###import matplotlib
###matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
###matplotlib.rcParams['xtick.direction'] = 'in'
###matplotlib.rcParams['ytick.direction'] = 'in'

# from sympy import *
###########Class and Function defined by HDPC################################
#sys.path.append("/sharefs/hbkg/user/nangyi/lick_tools")  # lib path
###sys.path.append("/sharefs/hbkg/user/luoqi/GRB/work/ihep4/lick_tools")
sys.path.append(os.path.abspath(os.path.dirname(__file__))+'/HXMT_GPS_py3Tools')
print('sucess2')
#import readxml
print('sucess3')
import LoadsrcLog
print('sucess3.0')
import funcLoadMod
print('sucess3.1')
from readxml import *
print('sucess3.2')
from LoadsrcLog import gensrctxt
import createpha_beta as createpha
print('sucess4')
import src_map as srcmap
import funchxmt_paraodic as hxmtpsf
from funcLoadMod import *
print('sucess5')
###import py3LoadMod
import xspec
from fitdata_plot import *
from mkdir_tree import *
print('sucess6')
print("import class success")


xspec.AllModels
def selfInitModel(src,outfilePri):
    print("New pot Runing >>>>>>>>>>> %s >>>>>>>>>>>\n"%src[0])
    time1= time.time()
    #from xspec import *
    #Model('grbm')
    xspec.AllModels
    #print(AllModels(1)(4))
    xspec.AllData.clear()
    xspec.AllModels.clear()
    xspec.AllModels.addPyMod(hxmtpsf.hxmtpsf, hxmtpsf.myModelParInfo, 'add')
    xspec.AllData("1:1 %s0.grp 2:2 %s1.grp 3:3 %s2.grp"%(outfilePri,outfilePri,outfilePri))
    d1 = xspec.AllData(1)
    d2 = xspec.AllData(2)
    d3 = xspec.AllData(3)
    xspec.AllModels += "powerlaw+hxmtpsf"
    m1 = xspec.AllModels(1)
    m2 = xspec.AllModels(2)
    m3 = xspec.AllModels(3)
    m1.setPars({1:"0 -1",2:"0. 0.1"})
    m2.setPars({1:"0 -1",2:"0. 0.1"})
    m3.setPars({1:"0 -1",2:"0. 0.1"})
    Fit.query='yes'
    ra,dec,norm = src[1],src[2],src[0]
    print(src,norm)
    m1.setPars({3:" 1          -1   "})
    m1.setPars({4:"60    -0.001"})
    m1.setPars({5:"%f,-0.01"%(ra)})
    m1.setPars({6:"%f,-0.01"%(dec)})
    m1.setPars({7:"%f,0.01"%(int(norm*18./100.) if int(norm*18./100.)<1000. else 900.)})
    m2.setPars({4:"0 -0.001"})
    m3.setPars({4:"-60  -0.001"})
    Fit.perform()
    Fit.error("max 100 1.0 1-7")
    tpchi = Fit.statistic
    tpdof = Fit.dof
    try:
        snr = 2*m1(7).values[0]/abs(m1(7).error[1]-m1(7).error[0])
        norm = m1(7).values[0]
    except:
        snr = -1
        norm = -1
    m1.setPars({7:"0.,-0.01"})
    zerochi = Fit.statistic
    deltachi = abs(tpchi - zerochi)
    low = m1(7).error[0]
    upp = m1(7).error[1]
    xspec.AllModels.clear()
    xspec.AllData.clear()
    print("New pot END!!!!!!!!!!!! %s >>>>>>>>>>>\n"%src[0])
    time2 = time.time()
    #del AllChains,AllData,AllModels,Fit,Plot,Xset
    return [ra,dec,norm,tpchi,tpdof,low,upp,deltachi]

def lcfit(cfg,program_tree,scan_tree,time1,Analogfile,MidductPath):
    with open("process", 'a') as f:
        f.write("Start\t---> " + str(time.strftime('%H:%M.%S, %Y.%m.%d', time.localtime(time.time()))) + "\n")
    #################################oooooo#####################################
    lcfile0, lcfile1, lcfile2, attfile, infilePri, outfilePri = GetXmlInf(cfg)

    print('infilePri:',infilePri)
    print('outfilePri:',outfilePri)

    createpha.cpha(lcfile0, "%s0" % outfilePri[:-3])
    createpha.cpha(lcfile1, "%s1" % outfilePri[:-3])
    createpha.cpha(lcfile2, "%s2" % outfilePri[:-3])
    time.sleep(2)

    AllData.clear()
    AllChains.clear()
    #del AllChains, AllData, AllModels, Fit, Plot, Xset
    #################################oooooo#####################################
    # crab_unit = 1000
    # bkg_flux = 10
    flux_thres = 0.001
    srcs, ralb, rahb, declb, dechb, ra_centre, dec_centre, srcsname, srcsnum, src_map_near, src_name_near = srcmap.src_new(
        cfg)
    print("srcs:::::::\n", srcs, '\n', ralb, rahb, declb, dechb, ra_centre, dec_centre)
    time.sleep(5)
    #################################oooooo#####################################
    sinfo = []
    for i, srcit in enumerate(srcs):
        info_tp = selfInitModel(srcit, outfilePri[:-3])
        sinfo.append(info_tp)
        if info_tp[5] < 0.001:
            srcs[i, 0] = 0.
        else:
            srcs[i, 0] = info_tp[2]
        if info_tp[3] == -1:
            srcs[i, 0] = -1
        if info_tp[6] > 20000:
            srcs[i, 0] = -1

    srcs_bright = srcs[srcs[:, 0] > flux_thres]
    srcsbrtname = srcsname[srcs[:, 0] > flux_thres]
    srcsbrtnum = srcsnum[srcs[:, 0] > flux_thres]

    srcs_zeros = srcs[srcs[:, 0] == 0.]
    srcszroname = srcsname[srcs[:, 0] == 0.]
    srcszronum = srcsnum[srcs[:, 0] == 0.]

    print("srcsbright::::", len(srcs_bright), '\n', srcsbrtname, srcs_bright, "\n")
    #from xspec import *

    ################################# Xset , loadmodel! ###################################
    try:
        logFile = Xset.openLog(outfilePri + '.log')
    except:
        print("Open log or par file failed!")
        sys.exit(1)

    xspec.AllModels.addPyMod(hxmtpsf.hxmtpsf, hxmtpsf.myModelParInfo, 'add')
    print("add hxmtpsf to the model")
    time.sleep(5)
    xspec.AllModels.clear()
    # Xset.save("hxmtmod.xcm",info='a')
    # Xset.parallel.leven = 14
    ###############################oooooooo#####################################
    temp = "powerlaw"
    k = 0
    for i, src in enumerate(srcs_bright):
        temp += "+hxmtpsf"
        k = k + 1
        print("srcs_bright number : ", k, temp)

    temp += "+hxmtpsf"

    parnumb = 5 * k + 2
    srclist = []
    goodra = []
    gooddec = []
    goodnorm = []
    goodsnr = []

    atthd = pf.open(attfile)
    atttb = atthd[1].data
    atttm = atttb.field(0)
    rapnt = float(atthd[0].header['RA_PNT'])
    decpnt = float(atthd[0].header['DEC_PNT'])
    dtobs = atthd[0].header['DATE-OBS']

    lchd = pf.open(lcfile2)
    lctb = lchd[1].data
    lctm = lctb.field(0)
    a = ((lctm[1:] - lctm[:-1]) != 1)
    c = np.nonzero(a)[0]
    pot = np.r_[0, c, lctm.shape[0] - 1]
    c = np.zeros(len(atttm))
    c = (c == 1)
    for i in range(len(pot) - 1):
        c = (c | ((atttm >= lctm[pot[i] + 1]) & (atttm <= lctm[pot[i + 1]])))

    tbatt = atttb[c]
    ratp = tbatt.field(1)
    dectp = tbatt.field(2)

    ###############load data and fit###########################################
    print("Fit Begin")
    AllData("1:1 %s0.grp 2:2 %s1.grp 3:3 %s2.grp" % (outfilePri[:-3], outfilePri[:-3], outfilePri[:-3]))
    time.sleep(5)
    d1 = AllData(1)
    d2 = AllData(2)
    d3 = AllData(3)
    xspec.AllModels += temp
    m1 = xspec.AllModels(1)
    m2 = xspec.AllModels(2)
    m3 = xspec.AllModels(3)
    m1.setPars({1: "0 -1", 2: "10 0.1"})
    m2.setPars({1: "0 -1", 2: "10 0.1"})
    m3.setPars({1: "0 -1", 2: "10 0.1"})
    for i, src in enumerate(srcs_bright):
        ra, dec, norm = src[1], src[2], src[0]
        print(ra, dec)
        m1.setPars({i * 5 + 3: " 1          -1   "})
        m1.setPars({i * 5 + 4: "60    -0.001"})
        m1.setPars({i * 5 + 5: "%f,-1.0" % ra})  # ,%f,%f,%f,%f"%(ra,Lrarg,Lrarg,Hrarg,Hrarg)})
        m1.setPars({i * 5 + 6: "%f,-1.0" % dec})  # ,%f,%f,%f,%f"%(dec,Ldecrg,Ldecrg,Hdecrg,Hdecrg)})
        m1.setPars({i * 5 + 7: "%f,0.5,%f,%f,%f,%f" % (norm if norm < 1000. else 900., 0., 0., 100000., 100000.)})
        m2.setPars({i * 5 + 4: "0 -0.001"})
        m3.setPars({i * 5 + 4: "-60  -0.001"})
        # m2(i*5+7).link = '1.0517 * %i'%(i*5+7)
        # m3(i*5+7).link = '1.0325 * %i'%(i*5+7)

    m1.setPars({k * 5 + 3: " 1          -1   "})
    m1.setPars({k * 5 + 4: "60    -0.001"})
    m1.setPars({k * 5 + 5: "0.,-1.0"})  # ,%f,%f,%f,%f"%(ra,Lrarg,Lrarg,Hrarg,Hrarg)})
    m1.setPars({k * 5 + 6: "0.,-1.0"})  # ,%f,%f,%f,%f"%(dec,Ldecrg,Ldecrg,Hdecrg,Hdecrg)})
    m1.setPars({k * 5 + 7: "%f,-0.01,%f,%f,%f,%f" % (0., -100000., -100000., 100000., 100000.)})
    m2.setPars({k * 5 + 4: "0 -0.001"})
    m3.setPars({k * 5 + 4: "-60  -0.001"})

    xspec.AllModels.show()
    Fit.query = 'yes'
    Fit.perform()
    # Fit.error("maximum 100 1.0 1-%s"%parnumb)
    # n1 = parnumb+7
    # n2 = 2*parnumb+12
    # Fit.error("maximum 100 1.0 %s"%n1)
    # Fit.error("maximum 100 1.0 %s"%n2)
    dataplot1(d1, d2, d3, m1, m2, m3, infilePri)

    Fit.error("maximum 100 1.0 1-%s" % parnumb)
    n1 = parnumb + 7
    n2 = 2 * parnumb + 12
    Fit.error("maximum 100 1.0 %s" % n1)
    Fit.error("maximum 100 1.0 %s" % n2)

    ########################################################################
    # Xset.parallel.reset()
    with open(outfilePri + '.src', "w") as srcFile:
        print("haha")
    with open(outfilePri + "New.src", "w") as newsrcFile:
        print("haha")
    ########################################################################
    print("****************** Start Error!*********************")
    with open("process", 'a') as f:
        f.write("Start Error  --->" + "\t" + str(time.strftime('%H:%M.%S, %Y.%m.%d', time.localtime(time.time()))) + "\n")
    red_chi = Fit.statistic / Fit.dof

    Ynew = []
    lenth = len(d1.noticed)
    lis = range(1, lenth, 100)
    for i in lis:
        AllData.ignore("**")
        AllData.notice("%s-%s" % (i, i + 100))
        Ynew.append(Fit.statistic)

    AllData.ignore("**")
    AllData.notice("%s-%s" % (lenth - 100, lenth))
    Ynew.append(Fit.statistic)

    AllData.notice("**")
    Ynew = (np.array(Ynew) > 500).sum()
    ########################################################################
    bkg1, bkg2, bkg3 = m1(2).values[0], m2(2).values[0], m3(2).values[0]
    upper1, upper2, upper3 = m1(2).error[1], m2(2).error[1], m3(2).error[1]
    low1, low2, low3 = m1(2).error[0], m2(2).error[0], m3(2).error[0]
    np.savetxt('bkg_Info.txt', [bkg1, bkg2, bkg3, low1, low2, low3, upper1, upper2, upper3])
    np.savetxt('all_psfInfo.txt', [lctm, m1.folded(1), m2.folded(2), m3.folded(3)])

    ########################################################################
    src_info = []
    src_cntr = []
    srcpot = []
    zerolist = []
    for i in range(7, parnumb + 1, 5):
        if m1(i).values[0] > .01:
            srclist.append(i)
        else:
            zerolist.append(i)

    for i in srclist:
        with open("process", 'a') as f:
            f.write("    --->" + str((i - 2) / 5) + "/" + str(len(srclist)) + "\t" + str(
                time.strftime('%H:%M.%S, %Y.%m.%d', time.localtime(time.time()))) + "\n")
        norm = m1(i).values[0]
        err1 = m1(i).error[0]
        err2 = m1(i).error[1]
        if err1 > err2:
            continue
        if err1 < 0.001:
            err1 = 2 * norm - err2
            snr = 2 * norm / abs(err2 - err1)
        else:
            snr = 2 * norm / abs(err2 - err1)

        if ((snr >= 5) & (snr < 9999999)):
            goodra.append(m1(i - 2).values[0])
            gooddec.append(m1(i - 1).values[0])
            goodnorm.append(m1(i).values[0])
            goodsnr.append(snr)

        tempmi = m1(i).values[0]
        tpchi = Fit.statistic
        fdtemp1 = np.array(m1.folded(1))
        fdtemp2 = np.array(m2.folded(2))
        fdtemp3 = np.array(m3.folded(3))
        ranow = m1(i - 2).values[0]
        decnow = m1(i - 1).values[0]
        m1.setPars({i: "0 -1"})
        newchi = Fit.statistic
        with open(outfilePri + '.src', "a") as srcFile:
            print(infilePri + "\t%s" % ((i - 2) / 5) + "\t" + srcsbrtname[int((i - 2) / 5 - 1)] + "\t" + str(
                srcsbrtnum[int((i - 2) / 5 - 1)] + 1) + "\n  Normal:\t%s" % norm + "\tLow:\t%s" % err1 + "\tUp:\t%s" % err2 + "\n  SNR:\t%s" % (
                snr) + "\t Chi-contribution:\t%s" % (tpchi - newchi) + "\n  Ra:\t%s" % (ranow) + "\n  Dec:\t%s" % (decnow),file=srcFile)
        src_info.append(
            [(i - 2) / 5, srcsbrtname[int((i - 2) / 5 - 1)], srcsbrtnum[int((i - 2) / 5 - 1)] + 1, ranow, decnow, norm,
             err1, err2])

        m1(i).frozen = True
        m1(i - 1).frozen = True
        m1(i - 2).frozen = True
        contri1 = fdtemp1 - np.array(m1.folded(1))
        contri2 = fdtemp2 - np.array(m2.folded(2))
        contri3 = fdtemp3 - np.array(m3.folded(3))
        src_cntr.append([contri1, contri2, contri3])
        # np.savetxt('src%s_psfInfo.txt'%(i-2)/5,[contri1,contri2,contri3])
        times1 = np.array(d1.noticed)[((fdtemp1 - np.array(m1.folded(1))) > 3)]
        times2 = np.array(d2.noticed)[((fdtemp2 - np.array(m2.folded(2))) > 3)]
        times3 = np.array(d3.noticed)[((fdtemp3 - np.array(m3.folded(3))) > 3)]
        m1.setPars({i: "%s -1" % tempmi})
        try:
            a = ((times1[1:] - times1[:-1]) != 1)
            c = np.nonzero(a)[0]
            pot1 = np.r_[0, c + 1]
            pot2 = np.r_[c, times1.shape[0] - 1]
            pot1 = times1[pot1]
            pot2 = times1[pot2]
            pot = (pot1 + pot2) / 2
            srcpot.append(pot)
        except:
            srcpot.append([])

        srcdtsum = 0
        srcmdsum = 0
        for j in range(len(pot1)):
            srcdtsum += (np.array(d1.values)[pot1[j]:pot2[j] + 1] - m1(2).values[0]).sum()
            srcmdsum += (contri1[pot1[j]:pot2[j] + 1]).sum()

        try:
            a = ((times2[1:] - times2[:-1]) != 1)
            c = np.nonzero(a)[0]
            pot1 = np.r_[0, c + 1]
            pot2 = np.r_[c, times2.shape[0] - 1]
            pot1 = times2[pot1]
            pot2 = times2[pot2]
            pot = (pot1 + pot2) / 2
            srcpot.append(pot)
        except:
            srcpot.append([])

        for j in range(len(pot1)):
            srcdtsum += (np.array(d2.values)[pot1[j]:pot2[j] + 1] - m2(2).values[0]).sum()
            srcmdsum += (contri2[pot1[j]:pot2[j] + 1]).sum()

        try:
            a = ((times3[1:] - times3[:-1]) != 1)
            c = np.nonzero(a)[0]
            pot1 = np.r_[0, c + 1]
            pot2 = np.r_[c, times3.shape[0] - 1]
            pot1 = times3[pot1]
            pot2 = times3[pot2]
            pot = (pot1 + pot2) / 2
            srcpot.append(pot)
        except:
            srcpot.append([])

        for j in range(len(pot1)):
            srcdtsum += (np.array(d3.values)[pot1[j]:pot2[j] + 1] - m3(2).values[0]).sum()
            srcmdsum += (contri3[pot1[j]:pot2[j] + 1]).sum()

        if i < 5 * k + 3:
            with open(outfilePri + '.src', "a") as srcFile:
                print("  Chi-Pct:\t%0.1f\t%0.1f" % (srcdtsum, srcmdsum * 100. / srcdtsum) + "%",file=srcFile)
        else:
            with open(outfilePri + "New.src", "a") as newsrcFile:
                print("  Chi-Pct:\t%0.1f\t%0.1f" % (srcdtsum, srcmdsum * 100. / srcdtsum) + "%",file=newsrcFile)

        if ((len(times1) > 0) | (len(times2) > 0) | (len(times3) > 0)):
            if i > 5 * k + 3:
                labelnm = int((i - 2) / 5)
            else:
                labelnm = srcsbrtname[int((i - 2) / 5) - 1]

            dataplot2(d1, d2, d3, m1, m2, m3, contri1, contri2, contri3, infilePri, i, labelnm)
        m1(i).frozen = True
        m1(i - 1).frozen = True
        m1(i - 2).frozen = True

    # Xset.closeLog()

    m1(2).frozen = True
    m2(2).frozen = True
    m3(2).frozen = True

    for i in zerolist:
        m1(i).frozen = True
        m1(i - 1).frozen = True
        m1(i - 2).frozen = True

    for i in zerolist:
        m1.setPars({i: "%f,0.5,%f,%f,%f,%f" % (0., -100000., -100000., 100000., 100000.)})
        try:
            Fit.perform()
        except:
            m1(i).frozen = True
            continue
        Fit.error('max 50. 1. %s' % i)
        norm = m1(i).values[0]
        err1 = m1(i).error[0]
        err2 = m1(i).error[1]
        if err1 < err2:
            snr = 2 * norm / abs(err2 - err1)
        else:
            snr = -1.
            m1(i).frozen = True
            continue
        ranow = m1(i - 2).values[0]
        decnow = m1(i - 1).values[0]
        with open(outfilePri + '.src', "a") as srcFile:
            print(infilePri + "\t%s" % ((i - 2) / 5) + "\t" + srcsbrtname[int((i - 2) / 5 - 1)] + "\t" + str(
                srcsbrtnum[int((
                                       i - 2) / 5 - 1)] + 1) + "\n  Normal:\t%s" % norm + "\tLow:\t%s" % err1 + "\tUp:\t%s" % err2 + "\n  SNR:\t%s" % (
                snr) + "\t Chi-contribution:\t%s" % (0.) + "\n  Ra:\t%s" % (ranow) + "\n  Dec:\t%s" % (decnow),file=srcFile)
        with open(outfilePri + '.src', "a") as srcFile:
            print("  Chi-Pct:\t%s\t%s" % ('Null', 'Null') + "%",file=srcFile)
        src_info.append(
            [(i - 2) / 5, srcsbrtname[int((i - 2) / 5 - 1)], srcsbrtnum[int((i - 2) / 5 - 1)] + 1, ranow, decnow, norm,
             err1, err2])
        fdtemp1 = np.array(m1.folded(1))
        fdtemp2 = np.array(m2.folded(2))
        fdtemp3 = np.array(m3.folded(3))
        m1.setPars({i: "%f,-0.5,%f,%f,%f,%f" % (0., -100000., -100000., 100000., 100000.)})
        contri1 = fdtemp1 - np.array(m1.folded(1))
        contri2 = fdtemp2 - np.array(m2.folded(2))
        contri3 = fdtemp3 - np.array(m3.folded(3))
        src_cntr.append([contri1, contri2, contri3])
        # np.savetxt('src%s_psfInfo.txt'%i,[contri1,contri2,contri3])

    for idx, src in enumerate(srcs_zeros):
        ranow, decnow, norm = src[1], src[2], src[0]
        i = k * 5 + 7
        m1.setPars({k * 5 + 5: "%f,-1.0" % ranow})  # ,%f,%f,%f,%f"%(ra,Lrarg,Lrarg,Hrarg,Hrarg)})
        m1.setPars({k * 5 + 6: "%f,-1.0" % decnow})  # ,%f,%f,%f,%f"%(dec,Ldecrg,Ldecrg,Hdecrg,Hdecrg)})
        m1.setPars({k * 5 + 7: "%f,0.01,%f,%f,%f,%f" % (0., -100000., -100000., 100000., 100000.)})
        try:
            Fit.perform()
        except:
            m1(i).frozen = True
            continue
        Fit.error('max 50. 1. %s' % i)
        norm = m1(i).values[0]
        err1 = m1(i).error[0]
        err2 = m1(i).error[1]
        if err1 < err2:
            snr = 2 * norm / abs(err2 - err1)
        else:
            snr = -1.
            m1(i).frozen = True
            continue
        with open(outfilePri + '.src', "a") as srcFile:
            print(infilePri + "\t%s" % (k + idx + 1) + "\t" + srcszroname[idx] + "\t" + str(srcszronum[
                                                                                                            idx] + 1) + "\n  Normal:\t%s" % norm + "\tLow:\t%s" % err1 + "\tUp:\t%s" % err2 + "\n  SNR:\t%s" % (
                snr) + "\t Chi-contribution:\t%s" % (0.) + "\n  Ra:\t%s" % (ranow) + "\n  Dec:\t%s" % (decnow),file=srcFile)
        with open(outfilePri + '.src', "a") as srcFile:
            print("  Chi-Pct:\t%s\t%s" % ('Null', 'Null') + "%",file=srcFile)
        src_info.append([idx + 1, srcszroname[idx], srcszronum[idx] + 1, ranow, decnow, norm, err1, err2])
        fdtemp1 = np.array(m1.folded(1))
        fdtemp2 = np.array(m2.folded(2))
        fdtemp3 = np.array(m3.folded(3))
        m1.setPars({i: "%f,-0.5,%f,%f,%f,%f" % (0., -100000., -100000., 100000., 100000.)})
        contri1 = fdtemp1 - np.array(m1.folded(1))
        contri2 = fdtemp2 - np.array(m2.folded(2))
        contri3 = fdtemp3 - np.array(m3.folded(3))
        src_cntr.append([contri1, contri2, contri3])

    np.save('src_Cntr.npy', src_cntr)
    #################################################################
    with open('src_Info.txt', 'w')as f:
        for ii in src_info:
            for ij in ii:
                f.write(str(ij) + '\t')
            f.write('\n')

    ###############################################################

    ratpanot = 0
    anotx = [0, 9999999999]
    anotxr = []
    anote = [0]
    anoter = []
    tplen = (len(srcpot) / 3)
    print('tplen:###########',tplen)
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
    dtemp_list = []
    for j in range(3):
        dtemp_list.append(xspec.AllModels(j + 1))

    dataplot3(d1, d2, d3, m1, m2, m3, infilePri, dtemp_list, srcpot, srclist)

    time2 = time.time()
    with open("process", 'a') as f:
        f.write("LcFit.png processed!\t---> %s" % ((time2 - time1) / 3600.0) + ' hours' + str(
            time.strftime('%H:%M.%S, %Y.%m.%d', time.localtime(time.time()))) + "\n")


    time2 = time.time()
    with open("process", 'a') as f:
        f.write("Test! 1\t---> %s" % ((time2 - time1) / 3600.0) + ' hours' + str(
            time.strftime('%H:%M.%S, %Y.%m.%d', time.localtime(time.time()))) + "\n")

    dataplot4(infilePri, srcs_bright, ratp, dectp, goodra, gooddec, goodsnr, goodnorm, src_map_near, time1)

    time2 = time.time()
    with open("process", 'a') as f:
        f.write("Test! 2\t---> %s" % ((time2 - time1) / 3600.0) + ' hours' + str(
            time.strftime('%H:%M.%S, %Y.%m.%d', time.localtime(time.time()))) + "\n")


    time2 = time.time()
    with open("process", 'a') as f:
        f.write("Test! 3\t---> %s" % ((time2 - time1) / 3600.0) + ' hours' + str(
            time.strftime('%H:%M.%S, %Y.%m.%d', time.localtime(time.time()))) + "\n")


    Xset.closeLog()
    atthd.close()
    lchd.close()

    time2 = time.time()
    print('fit finished ', time2)
    print("time expended: ", (time2 - time1) / 3600.0, ' hours')
    with open("process", 'a') as f:
        f.write("END\t---> %s" % ((time2 - time1) / 3600.0) + ' hours' + "\n")

    ########################gensrctxt########################
    gensrctxt(outfilePri[:-12], outfilePri[:-12])

    c = SkyCoord(ra=rapnt * u.degree, dec=decpnt * u.degree, frame='icrs')
    ldeg = float(c.galactic.l.degree)
    with open(Analogfile, 'a') as f:
        f.write(infilePri + '\t' + str(ldeg) + '\t' + str(dtobs) + '\t' + 'Y' + '\t' + 'Y\t' + str(
            time.strftime('%Y.%m.%d-%H:%M', time.localtime(time.time()))) + '\t' + str(Ynew) + '\n')

    os.system('cp src_P*.txt %s' % MidductPath)
    ##os.system('cp src_P*.txt %s' % MidductPath)
    print('Fit finished!')
    readcfg = loadDom(cfg)
    infilestr = readcfg.getTagText("infilenamelist")
    inpathstr = readcfg.getTagText("inpath")
    outpathstr = readcfg.getTagText("outpath")
    print(inpathstr, infilestr, outpathstr)
    inpathstr = inpathstr.strip()
    infilestr3 = infilestr.split()[3]
    infile3 = (inpathstr + infilestr3)  # .encode()
    ObsID = infile3.split("/")[-1][:13]
    fProd_obspath = scan_tree + '/HE/Prod/' + ObsID  ###
    mkdir_try(fProd_obspath)
    os.system('cp src_P*.txt %s/src_%s.txt'%(fProd_obspath,ObsID))
    os.system('python %s/HXMT_GPS_Program/Fit_Info.py %s'%(program_tree,cfg))


    with open(program_tree+'/gps_mission/lcfit_success.txt','a+') as f:
        print(outfilePri,file=f)
        print(infilePri,file=f)

def funcfit(cfg):
    #######read input output template paths and files from  config file##########
    dir_strlist = []
    try:
        with open('./dir_config.txt', 'r', encoding='utf-8') as f:
            for line in f.readlines():
                line = line.strip()
                dir_strlist.append(line)
        print(dir_strlist)
        program_tree = dir_strlist[0]
        scan_tree = dir_strlist[1]
    except:
        program_tree = '/sharefs/hbkg/user/luoqi/HXMT_SCAN/HE'
        scan_tree = '/sharefs/hbkg/data/SCAN'  # /sharefs/hbkg/user/luoqi/GRB/work/ihep4/SCAN'

    print('start main')

    time1 = time.time()
    print('cfg:', cfg)
    if len(cfg)>13:
        cfg_dir = cfg[:-14]
    else:
        cfg_dir = './'

    print('cfg_dir:',cfg_dir)
    print('change to cfg_dir==========>')
    os.chdir(cfg_dir)


    Analogfile = scan_tree + "/HE/Midd/Data_analys_log.txt"
    MidductPath = scan_tree + "/HE/Midd/srclist/"
    lcfit(cfg,program_tree,scan_tree,time1,Analogfile,MidductPath)
