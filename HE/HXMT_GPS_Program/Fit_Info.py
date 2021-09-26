from astropy.io import fits as pf
import sys,os,time
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
from mkdir_tree import *

###sys.path.append("/sharefs/hbkg/user/nangyi/lick_tools")#lib path
sys.path.append(os.path.abspath(os.path.dirname(__file__))+'/HXMT_GPS_py3Tools')
from readxml import *
print(sys.argv[0],sys.argv[1])
cfg = sys.argv[1]

def get_info(inattfile):
    print(inattfile)
    atthd = pf.open(inattfile)
    atttb = atthd[1].data
    atttm = atttb.field(0)
    rapnt = float(atthd[0].header['RA_PNT'])
    decpnt= float(atthd[0].header['DEC_PNT'])
    dtobs = atthd[0].header['DATE-OBS']
    c = SkyCoord(ra=rapnt*u.degree, dec=decpnt*u.degree, frame='icrs')
    ldeg = float(c.galactic.l.degree)
    bdeg = float(c.galactic.b.degree)
    t = Time(dtobs, format='isot', scale='utc')
    mjd = t.mjd
    atthd.close()
    del atthd,atttm
    return int(round(ldeg)),int(round(bdeg)),round(mjd,2),dtobs
########################################################################
readcfg = loadDom(cfg)
infilestr = readcfg.getTagText("infilenamelist")
inpathstr = readcfg.getTagText("inpath")
outpathstr = readcfg.getTagText("outpath")
outfilestr = readcfg.getTagText("outfilename")
print(inpathstr,infilestr,outpathstr)
inpathstr = inpathstr.strip()
infilestr0 = infilestr.split()[0]
infilestr1 = infilestr.split()[1]
infilestr2 = infilestr.split()[2]
infilestr3 = infilestr.split()[3]
outpathstr = outpathstr.split()[0]###
infile0 = (inpathstr+infilestr0)#.encode()
infile1 = (inpathstr+infilestr1)#.encode()
infile2 = (inpathstr+infilestr2)#.encode()
infile3 = (inpathstr+infilestr3)#.encode()
filename = infile3.split("/")[-1][:13]

os.chdir(outpathstr)
#------------------------------------------------------------
ldeg,bdeg,mjd,dtobs = get_info(infile3)
lctm,all1,all2,all3=np.loadtxt('all_psfInfo.txt')
fittm = time.ctime(os.stat("all_psfInfo.txt").st_mtime)
bkg1,bkg2,bkg3,low1,low2,low3,upper1,upper2,upper3 = np.loadtxt('bkg_Info.txt')
prihdr = pf.Header()
prihdr['TELESCOP'] = ('HXMT    ','Name of Telescope')
prihdr['INSTRUME'] = 'HE      '
prihdr['OBS_MODE'] = 'SAS    '
prihdr['OBS_ID'] = ('%s'%filename,'Proposition Number')
prihdr['OBS_MJD'] = mjd
prihdr['OBS_DATA'] = dtobs
prihdr['PNT_B'] = (ldeg,'[deg],Galactic Latitude')
prihdr['PNT_L'] = (bdeg,'[deg],Galactic Longitude')
prihdr['TSTART'] = (lctm[0],'Start Time')
prihdr['TSTOP'] = (lctm[-1],'Stop Time')
prihdr['SOFTWARE'] = ('hxmtsoftV2','light curve tools')
prihdr['TFITTED'] = (fittm,'Fiting date')

#prihdr['PSF'] = ('2018.01','PSF Calibrate Time')

prihdu = pf.PrimaryHDU(header=prihdr)
thdulist = pf.HDUList([prihdu])

numb,seq,names,ra,dec,norm,upper,lower = [float('nan'),float('nan'),float('nan')],[float('nan'),float('nan'),float('nan')],['BACKGRND_box1','BACKGRND_box2','BACKGRND_box3'],['nan','nan','nan'],['nan','nan','nan'],[bkg1,bkg2,bkg3],[upper1,upper2,upper3],[low1,low2,low3]
#upper[1:3],lower[1:3]= [float('nan'),float('nan')],[float('nan'),float('nan')]

with open('src_Info.txt','r')as f:
    for line in f:
        linelist = line.split('\t')
        numb.append(int(float(linelist[0])))###numb.append(int(linelist[0]))
        seq.append(int(float(linelist[2])))###seq.append(int(linelist[2]))
        names.append(linelist[1])
        ra.append(float(linelist[3]))
        dec.append(float(linelist[4]))
        norm.append(float(linelist[5]))
        upper.append(float(linelist[7]))
        lower.append(float(linelist[6]))

snr = 2.* np.array(norm)/(np.array(upper)- np.array(lower))
c1 = pf.Column(name='NAME',format='20A',array=names)
c2 = pf.Column(name='FITORDER',format='I',array=numb)
c3 = pf.Column(name='SEQUENSE',format='I',array=seq)
c4 = pf.Column(name='RA',format='E',array=ra)
c5 = pf.Column(name='DEC',format='E',array=dec)
c6 = pf.Column(name='NORM',format='E',array=norm)
c7 = pf.Column(name='NORM_UB',format='E',array=upper)
c8 = pf.Column(name='NORM_LB',format='E',array=lower)
c9 = pf.Column(name='SNR',format='E',array=snr)
coldefs = pf.ColDefs([c1, c2, c3, c4, c5, c6, c7, c8,c9])
tbhdu = pf.BinTableHDU.from_columns(coldefs,header = prihdr)
tbhdu.name = 'SOURCES'
thdulist.append(tbhdu)

lchd1 = pf.open(infile0)
lchd2 = pf.open(infile1)
lchd3 = pf.open(infile2)
tptime = lchd1[1].data.field(0)
yep = np.in1d(tptime,lctm)
lcdt1 = lchd1[1].data[yep]
lcdt2 = lchd2[1].data[yep]
lcdt3 = lchd3[1].data[yep]

src_cntr = np.load('src_Cntr.npy')
c1 = pf.Column(name='TIME',format='K',array=lctm)
c2 = pf.Column(name='ORG_CTS',format='E',array=lcdt1.field(1)+lcdt1.field(3))
c7 = pf.Column(name='SKY_BKG',format='E',array=lcdt1.field(3))
c3 = pf.Column(name='NET_CTS',format='E',array=lcdt1.field(1))
c4 = pf.Column(name='ORG_CTS_ERR',format='E',array=lcdt1.field(2))
c5 = pf.Column(name='PSF_SUM',format='E',array=all1)
c6 = pf.Column(name='RESIDUAL',format='E',array=lcdt1.field(1)-all1)
c8 = pf.Column(name='BKG_CONSTANT',format='E',array=np.array([norm[0]]*(len(lcdt1))))
coldefs = pf.ColDefs([c1, c2, c7, c3, c4, c5, c6,c8])
for i,it in enumerate(src_cntr):
    ct = pf.Column(name='%s'%names[3+i],format='E',array=src_cntr[i][0])
    coldefs.add_col(ct)

tbhdu = pf.BinTableHDU.from_columns(coldefs,header = prihdr)
tbhdu.name = 'DATA_BOX1'
thdulist.append(tbhdu)

c1 = pf.Column(name='TIME',format='K',array=lctm)
c2 = pf.Column(name='ORG_CTS',format='E',array=lcdt2.field(1)+lcdt2.field(3))
c7 = pf.Column(name='SKY_BKG',format='E',array=lcdt2.field(3))
c3 = pf.Column(name='NET_CTS',format='E',array=lcdt2.field(1))
c4 = pf.Column(name='ORG_CTS_ERR',format='E',array=lcdt2.field(2))
c5 = pf.Column(name='PSF_SUM',format='E',array=all2)
c6 = pf.Column(name='RESIDUAL',format='E',array=lcdt2.field(1)-all2)
c8 = pf.Column(name='BKG_CONSTANT',format='E',array=np.array([norm[1]]*(len(lcdt2))))
coldefs = pf.ColDefs([c1, c2, c7, c3, c4, c5, c6,c8])
for i,it in enumerate(src_cntr):
    ct = pf.Column(name='%s'%names[3+i],format='E',array=src_cntr[i][1])
    coldefs.add_col(ct)

tbhdu = pf.BinTableHDU.from_columns(coldefs,header = prihdr)
tbhdu.name = 'DATA_BOX2'
thdulist.append(tbhdu)

c1 = pf.Column(name='TIME',format='K',array=lctm)
c2 = pf.Column(name='ORG_CTS',format='E',array=lcdt3.field(1)+lcdt3.field(3))
c7 = pf.Column(name='SKY_BKG',format='E',array=lcdt3.field(3))
c3 = pf.Column(name='NET_CTS',format='E',array=lcdt3.field(1))
c4 = pf.Column(name='ORG_CTS_ERR',format='E',array=lcdt3.field(2))
c5 = pf.Column(name='PSF_SUM',format='E',array=all3)
c6 = pf.Column(name='RESIDUAL',format='E',array=lcdt3.field(1)-all3)
c8 = pf.Column(name='BKG_CONSTANT',format='E',array=np.array([norm[2]]*(len(lcdt3))))
coldefs = pf.ColDefs([c1, c2, c7, c3, c4, c5, c6,c8])
for i,it in enumerate(src_cntr):
    ct = pf.Column(name='%s'%names[3+i],format='E',array=src_cntr[i][2])
    coldefs.add_col(ct)

tbhdu = pf.BinTableHDU.from_columns(coldefs,header = prihdr)
tbhdu.name = 'DATA_BOX3'
thdulist.append(tbhdu)

atttm,attra,attdec,q1,q2,q3 = np.loadtxt('sate_Info.txt')
c1 = pf.Column(name='TIME',format='K',array=lctm)
c2 = pf.Column(name='RA',format='E',array=attra)
c3 = pf.Column(name='DEC',format='E',array=attdec)
c4 = pf.Column(name='Q1',format='E',array=q1)
c5 = pf.Column(name='Q2',format='E',array=q2)
c6 = pf.Column(name='Q3',format='E',array=q3)
coldefs = pf.ColDefs([c1, c2, c3, c4, c5, c6])
tbhdu = pf.BinTableHDU.from_columns(coldefs,header = prihdr)
tbhdu.name = 'ATT'
thdulist.append(tbhdu)

###bldhd = pf.open('/sharefs/hbkg/data/SCAN/nangyi/HE_data/lc_1811/lc/%sblind_g0_16.lc'%filename)
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
    scan_tree = '/sharefs/hbkg/data/SCAN'
bldhd = pf.open('%s/HE/Org/%s/%s_blind_g0_16.lc'%(scan_tree,filename,filename))
yep = np.in1d(bldhd[1].data.field(0),lctm)
blddt = bldhd[1].data[yep]
bldu = pf.BinTableHDU(data=blddt, header=prihdr, name= 'BlindDETECTOR')
thdulist.append(bldu)
thdulist.writeto('%s_FITDATA_HE.fits'%filename,overwrite = True)

fProd_obspath = scan_tree + '/HE/Prod/' + filename  ###
mkdir_try(fProd_obspath)
os.system('cp %s_FITDATA_HE.fits %s/%s_FITDATA_HE.fits'% (filename, fProd_obspath, filename))

'''
al = []
with open('log_info.txt','r')as f:
    for line in f:
        al.append(line[:8])

lt = []
with open('Data_analys_log.txt','r')as f:
    for line in f:
        nm = line[5:13]
        if nm not in al:
            commond = 'cd %s;pwd;cp /sharefs/hbkg/user/nangyi/code/Fit_Info.py ./;python Fit_Info.py config_he.xml;'%nm
            os.system(commond)
            al.append(nm)
            with open('log_info.txt','a')as fl:
                print>>fl,nm




'''
