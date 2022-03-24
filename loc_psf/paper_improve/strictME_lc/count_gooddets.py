import astropy.io.fits as pyfits
import numpy as np

c1_box0file = '/sharefs/hbkg/user/luoqi/psfl/genlc_strictME/7_40_1.0/lc/P021406400301me_g0_0-17.lc'
c2_box0file = '/sharefs/hbkg/user/luoqi/psfl/genlc_strictME/7_40_2.0/lc/P021406400301me_g0_0-17.lc'
c3_box0file = '/sharefs/hbkg/user/luoqi/psfl/genlc_strictME/7_40_3.0/lc/P021406400301me_g0_0-17.lc'
c4_box0file = '/sharefs/hbkg/user/luoqi/psfl/genlc_strictME/7_40_4.0/lc/P021406400301me_g0_0-17.lc'
c5_box0file = '/sharefs/hbkg/user/luoqi/psfl/genlc_strictME/7_40_5.0/lc/P021406400301me_g0_0-17.lc'
filelist = [c1_box0file,c2_box0file,c3_box0file,c4_box0file,c5_box0file]

def count_gdet(f):
    hdu = pyfits.open(f)
    pixel = hdu[-1].data['PIXEL']
    gdet = np.sum(pixel)
    return gdet

for i in range(5):
    f = filelist[i]
    gdet = count_gdet(f)
    print(i+1,gdet)