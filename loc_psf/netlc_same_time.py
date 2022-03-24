import numpy as np
from astropy.io import fits as pf
import sys,os
from glob import glob as glob
###sys.path.append("/hxmt/work/USERS/nangyi/tools")




heorigin_outpath = '/sharefs/hbkg/data/SCAN/luoqi/CrabLoc_HE/'
meorigin_outpath = '/sharefs/hbkg/data/SCAN/locME/'
leorigin_outpath = '/sharefs/hbkg/data/SCAN/luoqi/CrabLoc_LE/'
outpath= '/sharefs/hbkg/user/luoqi/psfl/crab_loc_nangyi/loc_luoqi'
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

heall0 = outpath+'/lqHEcrab_b0.fits'
heall1 = outpath+'/lqHEcrab_b1.fits'
heall2 = outpath+'/lqHEcrab_b2.fits'

meall0 = outpath+'/lqMEcrab_b0.fits'
meall1 = outpath+'/lqMEcrab_b1.fits'
meall2 = outpath+'/lqMEcrab_b2.fits'

leall0 = outpath +'/lqLE_SameTime_b0.fits'
#leall0 = outpath+'/lqLEcrab_b0.fits'
leall1 = outpath+'/lqLEcrab_b1.fits'
leall2 = outpath+'/lqLEcrab_b2.fits'

helist = [heall0,heall1,heall2]
melist = [meall0,meall1,meall2]
lelist = [leall0,leall1,leall2]

totlist = helist+melist+lelist

outHE0 = outpath +'/lqHE_SameTime_b0.fits'
outHE1 = outpath +'/lqHE_SameTime_b1.fits'
outHE2 = outpath +'/lqHE_SameTime_b2.fits'
outME0 = outpath +'/lqME_SameTime_b0.fits'
outME1 = outpath +'/lqME_SameTime_b1.fits'
outME2 = outpath +'/lqME_SameTime_b2.fits'
outLE0 = outpath +'/lqLE_SameTime_b0.fits'
outLE1 = outpath +'/lqLE_SameTime_b1.fits'
outLE2 = outpath +'/lqLE_SameTime_b2.fits'

outHElist = [outHE0,outHE1,outHE2]
outMElist = [outME0,outME1,outME2]
outLElist = [outLE0,outLE1,outLE2]

def he_write(tm,netct,neter,ct,cter,bkgct,bkger,outfits):
    col1 = pf.Column(name='Time', format='D', array=tm)
    col2 = pf.Column(name='Counts', format='D', array=netct)
    col3 = pf.Column(name='Stat_err', format='D', array=neter)
    col4 = pf.Column(name='Cts_org', format='D', array=ct)
    col5 = pf.Column(name='Cts_err', format='D', array=cter)
    col6 = pf.Column(name='Bkg', format='D', array=bkgct)
    col7 = pf.Column(name='Bkg_err', format='D', array=bkger)
    cols = pf.ColDefs([col1, col2, col3, col4, col5, col6, col7])
    tbhdu = pf.BinTableHDU.from_columns(cols)
    tbhdu.writeto(outfits, overwrite=True)

def me_write(tm,netct,neter,ct,bkgct,outfits):
    col1 = pf.Column(name='Time', format='D', array=tm)
    col2 = pf.Column(name='Counts', format='D', array=netct)
    col3 = pf.Column(name='Stat_err', format='D', array=neter)
    col4 = pf.Column(name='OGCounts', format='D', array=ct)
    col5 = pf.Column(name='BGCounts', format='D', array=bkgct)

    cols = pf.ColDefs([col1, col2, col3, col4, col5])
    tbhdu = pf.BinTableHDU.from_columns(cols)
    tbhdu.writeto(outfits, overwrite=True)

def hle_write(tm,netct,neter,ct,cter,bkgct,bkger,outfits):
    col1 = pf.Column(name='Time', format='D', array=tm)
    col2 = pf.Column(name='Counts', format='D', array=netct)
    col3 = pf.Column(name='Stat_err', format='D', array=neter)
    col4 = pf.Column(name='Cts_org', format='D', array=ct)
    col5 = pf.Column(name='Cts_err', format='D', array=cter)
    col6 = pf.Column(name='Bkg', format='D', array=bkgct)
    col7 = pf.Column(name='Bkg_err', format='D', array=bkger)
    cols = pf.ColDefs([col1, col2, col3, col4, col5, col6, col7])
    tbhdu = pf.BinTableHDU.from_columns(cols)
    tbhdu.writeto(outfits, overwrite=True)


def read_medata(file):
    hd = pf.open(file)
    tm = hd[1].data['TIME']
    netct = hd[1].data['Counts']
    neter = hd[1].data['Stat_err']
    ct = hd[1].data['OGCounts']
    bkgct = hd[1].data['BGCounts']
    return tm,netct,neter,ct,bkgct

def read_hledata(file):
    hd = pf.open(file)
    tm = hd[1].data['TIME']
    netct = hd[1].data['Counts']
    neter = hd[1].data['Stat_err']
    ct = hd[1].data['Cts_org']
    cter = hd[1].data['Cts_err']
    bkgct = hd[1].data['Bkg']
    bkger = hd[1].data['Bkg_err']
    return tm,netct,neter,ct,cter,bkgct,bkger



hetml = []
henetctl = []
heneterl = []
hectl = []
hecterl = []
hebkgctl = []
hebkgerl = []

for i in range(3):
    tm,netct,neter,ct,cter,bkgct,bkger = read_hledata(helist[i])
    hetml.append(tm)
    henetctl.append(netct)
    heneterl.append(neter)
    hectl.append(ct)
    hecterl.append(cter)
    hebkgctl.append(bkgct)
    hebkgerl.append(bkger)

metml = []
menetctl = []
meneterl = []
mectl = []
#hecterl = []
mebkgctl = []
#hebkgerl = []

for i in range(3):
    tm,netct,neter,ct,bkgct = read_medata(melist[i])
    metml.append(tm)
    menetctl.append(netct)
    meneterl.append(neter)
    mectl.append(ct)
    #hecterl.append(cter)
    mebkgctl.append(bkgct)
    #hebkgerl.append(bkger)



letml = []
lenetctl = []
leneterl = []
lectl = []
lecterl = []
lebkgctl = []
lebkgerl = []

for i in range(3):
    tm,netct,neter,ct,cter,bkgct,bkger = read_hledata(lelist[i])
    letml.append(tm)
    lenetctl.append(netct)
    leneterl.append(neter)
    lectl.append(ct)
    lecterl.append(cter)
    lebkgctl.append(bkgct)
    lebkgerl.append(bkger)


for j in range(0):

    maskhm = np.in1d(hetml[j],metml[j],invert=False)
    maskhl = np.in1d(hetml[j],letml[j],invert=False)
    maskall = np.logical_and(maskhl,maskhm)
    hetm = hetml[j][maskall]
    maskme = np.in1d(metml[j],hetm,invert=False)
    maskle0 = np.in1d(letml[0],hetm,invert=False)
    maskle1 = np.in1d(letml[1], hetm, invert=False)
    maskle2 = np.in1d(letml[2], hetm, invert=False)
    maskle = [maskle0,maskle1,maskle2]
    for i in range(3):
        try:
            hle_write(hetml[i][maskall], henetctl[i][maskall], heneterl[i][maskall], hectl[i][maskall], hecterl[i][maskall],
                      hebkgctl[i][maskall], hebkgerl[i][maskall], outHElist[i])
            me_write(metml[i][maskme], menetctl[i][maskme], meneterl[i][maskme], mectl[i][maskme],
                     mebkgctl[i][maskme], outMElist[i])
            # hle_write(letml[i][maskle], lenetctl[i][maskle], leneterl[i][maskle], lectl[i][maskle], lecterl[i][maskle],
            #           lebkgctl[i][maskle], lebkgerl[i][maskle], outLElist[i])
            hle_write(letml[i][maskle[i]], lenetctl[i][maskle[i]], leneterl[i][maskle[i]], lectl[i][maskle[i]], lecterl[i][maskle[i]],
                      lebkgctl[i][maskle[i]], lebkgerl[i][maskle[i]], outLElist[i])
        except:
            print(j,i)




