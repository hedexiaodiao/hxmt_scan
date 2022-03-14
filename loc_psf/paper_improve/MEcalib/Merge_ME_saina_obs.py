import numpy as np
from astropy.io import fits as pf
import sys,os
from glob import glob as glob
###sys.path.append("/hxmt/work/USERS/nangyi/tools")
from Fitstogether2 import fitstogetherMulti

listall = []
with open('041crab.all','r')as f:
    for i in f:
        listall.append(i[:13])


outpath= './'###'/sharefs/hbkg/data/SCAN/PSF/HE/all_2020/'

attall = []
lc0all = []
lc1all = []
lc2all = []

for i in range(len(listall)):
    ObsID = listall[i]
    stObsID = ObsID[0:8]

    attpath = '/sharefs/hbkg/user/saina/data294/P0101294/'
    lcpath = '/sharefs/hbkg/user/saina/data294/P0101294/'
    Wcutattfile = attpath+ObsID+'/ME/att.fits'

    attfile = Wcutattfile
    lc0file = lcpath+ObsID+'/ME/me_lc_box0_small_cut.fits'
    lc1file = lcpath+ObsID+'/ME/me_lc_box1_small_cut.fits'
    lc2file = lcpath+ObsID+'/ME/me_lc_box2_small_cut.fits'
    attflag = os.path.exists(attfile)
    lc0flag = os.path.exists(lc0file)
    lc1flag = os.path.exists(lc1file)
    lc2flag = os.path.exists(lc2file)
    if attflag and lc0flag and lc1flag and lc2flag:
        attall.append(attfile)
        lc0all.append(lc0file)
        lc1all.append(lc1file)
        lc2all.append(lc2file)
    else:
        print("lc file no exist:",ObsID,attflag,lc0flag,lc1flag,lc2flag)


fitstogetherMulti(lc0all,outpath+'me_b0.fits')
fitstogetherMulti(lc1all,outpath+'me_b1.fits')
fitstogetherMulti(lc2all,outpath+'me_b2.fits')
fitstogetherMulti(attall,outpath+'me_att.fits')


