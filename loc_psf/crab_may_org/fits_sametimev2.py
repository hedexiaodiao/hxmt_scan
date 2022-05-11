import numpy as np
from astropy.io import fits as pf

inpath= '/sharefs/hbkg/user/luoqi/psfl/crab_loc_v1/data_crab_v2'
outpath= '/sharefs/hbkg/user/luoqi/psfl/crab_loc_v1/data_crab_v2/sametime'

he_inlist = [inpath+'/P0101299allHE_b0.fits',inpath+'/P0101299allHE_b1.fits',inpath+'/P0101299allHE_b2.fits']
me_inlist = [inpath+'/P0101299allME_b0.fits',inpath+'/P0101299allME_b1.fits',inpath+'/P0101299allME_b2.fits']
le_inlist = [inpath+'/P0101299allLE_b0.fits',inpath+'/P0101299allLE_b1.fits',inpath+'/P0101299allLE_b2.fits']

he_outlist = [outpath+'/P0101299allHE_b0.fits',outpath+'/P0101299allHE_b1.fits',outpath+'/P0101299allHE_b2.fits']
me_outlist = [outpath+'/P0101299allME_b0.fits',outpath+'/P0101299allME_b1.fits',outpath+'/P0101299allME_b2.fits']
le_outlist = [outpath+'/P0101299allLE_b0.fits',outpath+'/P0101299allLE_b1.fits',outpath+'/P0101299allLE_b2.fits']

inlist = [he_inlist,me_inlist,le_inlist]
outlist = [he_outlist,me_outlist,le_outlist]

he0 = pf.open(he_inlist[0])
he1 = pf.open(he_inlist[1])
he2 = pf.open(he_inlist[2])
me0 = pf.open(me_inlist[0])
me1 = pf.open(me_inlist[1])
me2 = pf.open(me_inlist[2])
le0 = pf.open(le_inlist[0])
le1 = pf.open(le_inlist[1])
le2 = pf.open(le_inlist[2])
hepf_list = [he0,he1,he2]
mepf_list = [me0,me1,me2]
lepf_list = [le0,le1,le2]
pf_list = [hepf_list,mepf_list,lepf_list]


timehe0 = he0[1].data['Time']
timehe1 = he1[1].data['Time']
timehe2 = he2[1].data['Time']
timeme0 = me0[1].data['Time']
timeme1 = me1[1].data['Time']
timeme2 = me2[1].data['Time']
timele0 = le0[1].data['Time']
timele1 = le1[1].data['Time']
timele2 = le2[1].data['Time']

timehe_list = [timehe0,timehe1,timehe2]
timeme_list = [timeme0,timeme1,timeme2]
timele_list = [timele0,timele1,timele2]


for i in [0,1,2]:
    for j in range(3):
        k = j
        mask0 = np.in1d(timehe_list[k], timeme_list[k][np.in1d(timeme_list[k], timele_list[k])])
        mask1 = np.in1d(timeme_list[k], timele_list[k][np.in1d(timele_list[k], timehe_list[k])])
        mask2 = np.in1d(timele_list[k], timehe_list[k][np.in1d(timehe_list[k], timeme_list[k])])
        mask_list = [mask0, mask1, mask2]
        hdu = pf_list[i][j]
        tm = hdu[1].data['Time'][mask_list[i]]
        cts = hdu[1].data['Counts'][mask_list[i]]
        stat_err = hdu[1].data['Stat_err'][mask_list[i]]
        cts_org = hdu[1].data['Cts_org'][mask_list[i]]
        cts_err = hdu[1].data['Cts_err'][mask_list[i]]
        bkg = hdu[1].data['Bkg'][mask_list[i]]
        bkg_err = hdu[1].data['Bkg_err'][mask_list[i]]
        col1 = pf.Column(name='Time', format='D', array=tm)
        col2 = pf.Column(name='Counts', format='D', array=cts)
        col3 = pf.Column(name='Stat_err', format='D', array=stat_err)
        col4 = pf.Column(name='Cts_org', format='D', array=cts_org)
        col5 = pf.Column(name='Cts_err', format='D', array=cts_err)
        col6 = pf.Column(name='Bkg', format='D', array=bkg)
        col7 = pf.Column(name='Bkg_err', format='D', array=bkg_err)
        cols = pf.ColDefs([col1, col2, col3, col4, col5, col6, col7])
        tbhdu = pf.BinTableHDU.from_columns(cols)
        tbhdu.writeto(outlist[i][j])

# for i in [1]:
#     for j in range(3):
#         k = j
#         mask0 = np.in1d(timehe_list[k], timeme_list[k][np.in1d(timeme_list[k], timele_list[k])])
#         mask1 = np.in1d(timeme_list[k], timele_list[k][np.in1d(timele_list[k], timehe_list[k])])
#         mask2 = np.in1d(timele_list[k], timehe_list[k][np.in1d(timehe_list[k], timeme_list[k])])
#         mask_list = [mask0, mask1, mask2]
#         hdu = pf_list[i][j]
#         tm = hdu[1].data['Time'][mask_list[i]]
#         cts = hdu[1].data['Counts'][mask_list[i]]
#         stat_err = hdu[1].data['Stat_err'][mask_list[i]]
#         cts_org = hdu[1].data['OGCounts'][mask_list[i]]
#         bkg_org = hdu[1].data['BGCounts'][mask_list[i]]
#         col1 = pf.Column(name='Time', format='D', array=tm)
#         col2 = pf.Column(name='Counts', format='D', array=cts)
#         col3 = pf.Column(name='Stat_err', format='D', array=stat_err)
#         col4 = pf.Column(name='OGCounts', format='D', array=cts_org)
#         col5 = pf.Column(name='BGCounts', format='D', array=bkg_org)
#         cols = pf.ColDefs([col1, col2, col3, col4, col5])
#         tbhdu = pf.BinTableHDU.from_columns(cols)
#         tbhdu.writeto(outlist[i][j])
