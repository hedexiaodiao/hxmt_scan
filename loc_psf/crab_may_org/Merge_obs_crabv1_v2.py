import numpy as np
from astropy.io import fits as pf
import sys,os
import glob
###sys.path.append("/hxmt/work/USERS/nangyi/tools")
from Fitstogether2 import fitstogetherMulti

listall = []

id_end2 = [1,2,3,4,6,7,8,9,10,11,12,13]
for i in id_end2:
    for j in range(24):
        obsid = 'P0101299{:03d}{:02d}'.format(i,j)
        basedatapath = '/hxmt/work/HXMT-DATA/1L/'+"A%s/%s/%s/"%(obsid[1:3],obsid[:8],obsid[:11])
        ###if os.path.exists(basedatapath):
        ###    print(obsid,'exist!!!!')
        datapath_list = glob.glob(basedatapath+obsid+'-*')
        if len(datapath_list)>0:
            datapath = datapath_list[0]
            if os.path.exists(datapath):
                if i==8 and j==18: continue;
                if i==8 and j==19: continue;
                print(obsid,'exist!!!!')
                listall.append(obsid)

heorigin_outpath = '/sharefs/hbkg/data/SCAN/luoqi/Crabv1_HE/'
meorigin_outpath = '/sharefs/hbkg/data/SCAN/locv1/'
leorigin_outpath = '/sharefs/hbkg/data/SCAN/luoqi/CrabLoc_LE/'
outpath= '/sharefs/hbkg/user/luoqi/psfl/crab_loc_v1/data_crab_v1'
###------------
'''
att:
/sharefs/hbkg/data/SCAN/luoqi/CrabLoc_HE/att/P0101299008*crtatt

he:
/sharefs/hbkg/data/SCAN/luoqi/CrabLoc_HE/lc/P0101299008*hebox*.fits

me: 
/sharefs/hbkg/data/SCAN/locME/lc/P0101299008*mebox*.fits

le: 
/sharefs/hbkg/data/SCAN/luoqi/CrabLoc_LE/lc/P0101299008*lebox*.fits
'''


print(listall)

attall = []
helc0all = []
helc1all = []
helc2all = []

melc0all = []
melc1all = []
melc2all = []

lelc0all = []
lelc1all = []
lelc2all = []

for i in range(len(listall)):
    ObsID = listall[i]
    stObsID = ObsID[0:8]
    attfile = '/sharefs/hbkg/data/SCAN/luoqi/CrabLoc_HE/att/%scrtatt'%(ObsID)
    attflag = os.path.exists(attfile)

    helcfile = heorigin_outpath + 'lc/%she' % ObsID
    helc0file = helcfile+'box0.fits'
    helc1file = helcfile+'box1.fits'
    helc2file = helcfile+'box2.fits'
    helc0flag = os.path.exists(helc0file)
    helc1flag = os.path.exists(helc1file)
    helc2flag = os.path.exists(helc2file)

    melcfile = meorigin_outpath + 'lc/%sme' % ObsID
    melc0file = melcfile + 'box0.fits'
    melc1file = melcfile + 'box1.fits'
    melc2file = melcfile + 'box2.fits'
    melc0flag = os.path.exists(melc0file)
    melc1flag = os.path.exists(melc1file)
    melc2flag = os.path.exists(melc2file)

    lelcfile = leorigin_outpath + 'lc/%sle' % ObsID
    lelc0file = lelcfile + 'box0.fits'
    lelc1file = lelcfile + 'box1.fits'
    lelc2file = lelcfile + 'box2.fits'
    lelc0flag = os.path.exists(lelc0file)
    lelc1flag = os.path.exists(lelc1file)
    lelc2flag = os.path.exists(lelc2file)

    if attflag and helc0flag and helc1flag and helc2flag and melc0flag and melc1flag and melc2flag and lelc0flag and lelc1flag and lelc2flag:
        attall.append(attfile)
        helc0all.append(helc0file)
        helc1all.append(helc1file)
        helc2all.append(helc2file)

        melc0all.append(melc0file)
        melc1all.append(melc1file)
        melc2all.append(melc2file)

        lelc0all.append(lelc0file)
        lelc1all.append(lelc1file)
        lelc2all.append(lelc2file)
    else:
        print("lc file no exist:",ObsID,attflag,helc0flag,helc1flag,helc2flag,melc0flag,melc1flag,melc2flag,lelc0flag,lelc1flag, lelc2flag)

###print(attall)
###fitstogetherMulti(attall,outpath+'/P0101299all_att.fits')

print(helc0all,helc1all,helc2all)
fitstogetherMulti(helc0all,outpath+'/P0101299allHE_b0.fits')
fitstogetherMulti(helc1all,outpath+'/P0101299allHE_b1.fits')
fitstogetherMulti(helc2all,outpath+'/P0101299allHE_b2.fits')
print('he finish')

print(melc0all,melc1all,melc2all)
fitstogetherMulti(melc0all,outpath+'/P0101299allME_b0.fits')
fitstogetherMulti(melc1all,outpath+'/P0101299allME_b1.fits')
fitstogetherMulti(melc2all,outpath+'/P0101299allME_b2.fits')
print('me finish')

print(lelc0all,lelc1all,lelc2all)
fitstogetherMulti(lelc0all,outpath+'/P0101299allLE_b0.fits')
fitstogetherMulti(lelc1all,outpath+'/P0101299allLE_b1.fits')
fitstogetherMulti(lelc2all,outpath+'/P0101299allLE_b2.fits')
print('le finish')


