import numpy as np
from astropy.io import fits as pf
import sys,os
from glob import glob as glob
###sys.path.append("/hxmt/work/USERS/nangyi/tools")
from Fitstogether2 import fitstogetherMulti

listall = []
with open('crablist.txt','r')as f:
    for i in f:
        listall.append('P0'+i[:9]+'01')

outpath= '/sharefs/hbkg/user/luoqi/psfl/ME'

print(listall)

attall = []
lc0all = []
lc1all = []
lc2all = []
for i in range(len(listall)):
    ObsID = listall[i]
    stObsID = ObsID[0:8]
    attfile = '/sharefs/hbkg/user/saina/data294/%s/%s/ME/att.fits'%(stObsID,ObsID)
    lc0file = '/sharefs/hbkg/user/saina/data294/%s/%s/ME/me_lc_box0_small_cut.fits'%(stObsID,ObsID)
    lc1file = '/sharefs/hbkg/user/saina/data294/%s/%s/ME/me_lc_box1_small_cut.fits'%(stObsID,ObsID)
    lc2file = '/sharefs/hbkg/user/saina/data294/%s/%s/ME/me_lc_box2_small_cut.fits'%(stObsID,ObsID)
    if os.path.exists(attfile):
        attall.append(attfile)
    if os.path.exists(lc0file):
        lc0all.append(lc0file)
    if os.path.exists(lc1file):
        lc1all.append(lc1file)
    if os.path.exists(lc2file):
        lc2all.append(lc2file)

print(attall)
fitstogetherMulti(lc0all,outpath+'/me_b0.fits')
fitstogetherMulti(lc1all,outpath+'/me_b1.fits')
fitstogetherMulti(lc2all,outpath+'/me_b2.fits')
fitstogetherMulti(attall,outpath+'/me_att.fits')


