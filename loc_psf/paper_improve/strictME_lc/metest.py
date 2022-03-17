import sys
import numpy as np
import xlrd
import os
import glob
from astropy.io import fits as pf
#bdfile = '/home/hxmt/hxmtsoft2/hxmtsoftv2.02/install/x86_64-pc-linux-gnu-libc2.12/refdata/medetectorstatus.fits'
from tool import *
sys.path.append("/sharefs/hbkg/user/nangyi/HXMT_GPS_Tools")
from Lctools_202004 import correct

obsid = sys.argv[1]
pfile = "/home/hxmt/nangyi/pfiles/me%s"%str(obsid)
os.system('mkdir %s'%pfile);os.system("rm %s/*"%pfile);
pfilepath = "%s;/home/hxmt/hxmtsoft2/hxmtsoftv2.02/install/x86_64-pc-linux-gnu-libc2.12/syspfiles"%pfile
cmdinit = 'sleep 1;source /sharefs/hbkg/user/nangyi/hxmtsoft_v2.02.sh ;export PFILES="%s" '%(pfilepath)
minpi = 68
#maxpi = 166
#inpi = 166
maxpi = 631
basedatapath = '/hxmt/work/HXMT-DATA/1L/'+"A%s/%s/%s/"%(obsid[1:3],obsid[:8],obsid[:11])
datapath = glob.glob(basedatapath+obsid+'-*')[0]
basefile = get_rawdata(datapath, instrument="ME")

outpath = '/sharefs/hbkg/data/SCAN/nangyi/ME_data/'
pifile = outpath + 'pi/%smepi.fits'%obsid
grdfile = outpath + 'grd/%smegrd.fits'%obsid
deadfile = outpath + 'grd/%smedead.fits'%obsid
gtifile = outpath + 'gti/%smegti.fits'%obsid
newgtifile = outpath + 'gti/%smegtinew.fits'%obsid
newbdfile = outpath + 'gti/%smebd.fits'%obsid
screenfile = outpath + 'screen/%smescreen.fits'%obsid
lcfile = outpath + 'lc/%sme'%obsid
lcnametxt = outpath + 'lc/%smelc.txt'%obsid
bkgfile = outpath + 'lc/%smebkg'%obsid
#basefile['EHK']='/sharefs/hbkg/data/SCAN/nangyi/HE_point/ehk/%sehk.fits'%obsid

crtattfile = outpath + 'att/%scrtatt'%obsid
correct(basefile['ATT'],crtattfile)

cmd=''
if not os.path.exists(pifile):
    cmd='; mepical evtfile=%s tempfile=%s outfile=%s'%(basefile['EVT'],basefile['TH'],pifile)

if not os.path.exists(grdfile):
    cmd +=' ;megrade evtfile=%s outfile=%s deadfile=%s binsize=1'%(pifile,grdfile,deadfile)

if not os.path.exists(gtifile):
    cmd +=' ;megtigen tempfile="%s" ehkfile="%s" outfile="%s" defaultexpr=NONE expr="SAA_FLAG==0&&T_SAA>100&&TN_SAA>100&&ANG_DIST<=360" clobber=yes history=yes'%(basefile['TH'], basefile['EHK'], gtifile)
    #cmd +=' ;megti %s %s %s $HEADAS/refdata/medetectorstatus.fits %s'%(grdfile, gtifile, newgtifile, newbdfile)

newbdfile = outpath + 'gti/mebadall.fits'

if not os.path.exists(screenfile):
    cmd +=' ;mescreen evtfile=%s gtifile=%s outfile=%s baddetfile=%s userdetid="0-53"'%(grdfile,gtifile,screenfile,newbdfile)

os.system(cmdinit + cmd)
thd = pf.open(screenfile)#,endcard=False)
try:
    print len(thd)
except:
    pass

startm = int(thd[2].data.field(0)[0])+1
stoptm = int(thd[2].data.field(1)[-1])-1

thd.close();del thd
#float('l')
cmd = '; melcgen evtfile=%s deadfile=%s outfile=%s userdetid="0-7,11-17;18-25,29-35;36-43,47-53" starttime=%s stoptime=%s minPI=%s maxPI=%s binsize=0.5 deadcorr=yes'%(screenfile,deadfile,lcfile,startm,stoptm,minpi,maxpi)

os.system(cmdinit + cmd)


#float('l')
cmd = ' ; ls -t %s_g0*.lc | head -1 >  %s'%(lcfile,lcnametxt)
cmd += ' ; python mebkgmap_2012005_lc_err.py lc %s %s %s %s %s %s %s %s %s %s'%(screenfile,basefile['EHK'],gtifile,deadfile,basefile['TH'],lcnametxt,minpi,maxpi,bkgfile,newbdfile)
os.system(cmdinit + cmd)

hd=pf.open(lcfile+'bkg_box1.lc')
bgtm = hd[1].data.field(0)
bgct = hd[1].data.field(1)
hd.close();del hd
a = peakdel(bgtm,bgct,50,15)
print 'haha'
gti = bgtm[a]
print 'haha'

gtifile_short = gtifile = outpath + 'gti/%smegtishort.fits'%obsid
screenfile_short = outpath + 'screen/%smescreenshort.fits'%obsid
cmd =' ;megtigen tempfile="%s" ehkfile="%s" outfile="%s" defaultexpr=NONE expr="ELV>8&&COR>8&&SAA_FLAG==0&&T_SAA>100&&TN_SAA>100&&ANG_DIST<=0.02" clobber=yes history=yes'%(basefile['TH'], basefile['EHK'], gtifile_short)
cmd +=' ;mescreen evtfile=%s gtifile=%s outfile=%s baddetfile=%s userdetid="0-53"'%(grdfile,gtifile_short,screenfile_short,newbdfile)
os.system(cmdinit + cmd)
thd = pf.open(screenfile)#,endcard=False)
try:
    print len(thd)
except:
    pass
startm = int(thd[2].data.field(0)[0])+1
stoptm = int(thd[2].data.field(1)[-1])-1
thd.close();del thd
cmd = '; melcgen evtfile=%s deadfile=%s outfile=%s userdetid="0-7,11-17;18-25,29-35;36-43,47-53" starttime=%s stoptime=%s minPI=%s maxPI=%s binsize=0.5 deadcorr=yes'%(screenfile_short,deadfile,lcfile+'gti',startm,stoptm,minpi,maxpi)
os.system(cmdinit + cmd)
hd=pf.open('%sgti_g0_0-17.lc'%lcfile)
gtigt = hd[1].data.field(0)
hd.close();del hd
#float('l')

sfile=pf.open("%s"%newbdfile)
bdnum=sfile[1].data.field(0)
bdbox0 = np.in1d(range(256)+range(352,576),bdnum).sum()
bdbox1 = np.in1d(range(576,832)+range(928,1152),bdnum).sum()
bdbox2 = np.in1d(range(1152,1408)+range(1504,1728),bdnum).sum()
sfile.close();del sfile
r0=480./(480-bdbox0)
r1=480./(480-bdbox1)
r2=480./(480-bdbox2)
rall=[r0,r1,r2]
#print np.array(range(256)+range(352,576))[~np.in1d(range(256)+range(352,576),bdnum)],np.array(range(256)+range(352,576))[~np.in1d(range(256)+range(352,576),bdnum)].shape
#r0,r1,r2=1,1,1

hd1=pf.open('%s_g0_0-17.lc'%lcfile)
hd2=pf.open('%s_box0.lc'%bkgfile)
yp = (hd1[1].data['FRACEXP']==1.) & np.in1d(hd1[1].data['TIME'],gti) & np.in1d(hd1[1].data['TIME'],gtigt)
tm = hd1[1].data['TIME'][yp]
ct = hd1[1].data['COUNTS'][yp]*r0
cter = hd1[1].data['ERROR'][yp]*r0
bkgct = hd2[1].data['RATE'][yp]*r0/2.
bkger = hd2[1].data['Error'][yp]*r0/2.
#print ct.sum()/bkgct.sum()
#with open('bkg0.txt','a')as f:
#    print>>f,ct.sum(),bkgct.sum(),obsid

netct = ct-bkgct
neter = np.sqrt(cter**2+bkger**2)
col1 = pf.Column(name = 'Time', format='D',array = tm)
col2 = pf.Column(name = 'Counts', format='D',array = netct)
col3 = pf.Column(name = 'Stat_err', format='D',array = neter)
#col4 = pf.Column(name = 'Cts_org', format='D',array = ct)
#col5 = pf.Column(name = 'Cts_err', format='D',array = cter)
#col6 = pf.Column(name = 'Bkg', format='D',array = bkgct)
#col7 = pf.Column(name = 'Bkg_err', format='D',array = bkger)
col4 = pf.Column(name = 'OGCounts', format='D',array = ct)
col5 = pf.Column(name = 'BGCounts', format='D',array = bkgct)

cols = pf.ColDefs([col1,col2,col3,col4,col5])
tbhdu = pf.BinTableHDU.from_columns(cols)
tbhdu.writeto(lcfile + "box0.fits",overwrite=True)
hd1.close();hd2.close();del hd1,hd2


hd1=pf.open('%s_g1_18-35.lc'%lcfile)
hd2=pf.open('%s_box1.lc'%bkgfile)
yp = (hd1[1].data['FRACEXP']==1.) & np.in1d(hd1[1].data['TIME'],gti) & np.in1d(hd1[1].data['TIME'],gtigt)
tm = hd1[1].data['TIME'][yp]
ct = hd1[1].data['COUNTS'][yp]*r1
cter = hd1[1].data['ERROR'][yp]*r1
bkgct = hd2[1].data['RATE'][yp]*r1/2.
bkger = hd2[1].data['Error'][yp]*r1/2.
netct = ct-bkgct
#print ct.sum()/bkgct.sum()
#with open('bkg1.txt','a')as f:
#    print>>f,ct.sum(),bkgct.sum(),obsid


neter = np.sqrt(cter**2+bkger**2)
col1 = pf.Column(name = 'Time', format='D',array = tm)
col2 = pf.Column(name = 'Counts', format='D',array = netct)
col3 = pf.Column(name = 'Stat_err', format='D',array = neter)
#col4 = pf.Column(name = 'Cts_org', format='D',array = ct)
#col5 = pf.Column(name = 'Cts_err', format='D',array = cter)
#col6 = pf.Column(name = 'Bkg', format='D',array = bkgct)
#col7 = pf.Column(name = 'Bkg_err', format='D',array = bkger)
col4 = pf.Column(name = 'OGCounts', format='D',array = ct)
col5 = pf.Column(name = 'BGCounts', format='D',array = bkgct)

cols = pf.ColDefs([col1,col2,col3,col4,col5])
tbhdu = pf.BinTableHDU.from_columns(cols)
tbhdu.writeto(lcfile + "box1.fits",overwrite=True)
hd1.close();hd2.close();del hd1,hd2


hd1=pf.open('%s_g2_36-53.lc'%lcfile)
hd2=pf.open('%s_box2.lc'%bkgfile)
yp = (hd1[1].data['FRACEXP']==1.) & np.in1d(hd1[1].data['TIME'],gti) & np.in1d(hd1[1].data['TIME'],gtigt)
tm = hd1[1].data['TIME'][yp]
ct = hd1[1].data['COUNTS'][yp]*r2
cter = hd1[1].data['ERROR'][yp]*r2
bkgct = hd2[1].data['RATE'][yp]*r2/2.
bkger = hd2[1].data['Error'][yp]*r2/2.
netct = ct-bkgct
#print ct.sum()/bkgct.sum()
#with open('bkg2.txt','a')as f:
#    print>>f,ct.sum(),bkgct.sum(),obsid


neter = np.sqrt(cter**2+bkger**2)
col1 = pf.Column(name = 'Time', format='D',array = tm)
col2 = pf.Column(name = 'Counts', format='D',array = netct)
col3 = pf.Column(name = 'Stat_err', format='D',array = neter)
#col4 = pf.Column(name = 'Cts_org', format='D',array = ct)
#col5 = pf.Column(name = 'Cts_err', format='D',array = cter)
#col6 = pf.Column(name = 'Bkg', format='D',array = bkgct)
#col7 = pf.Column(name = 'Bkg_err', format='D',array = bkger)
col4 = pf.Column(name = 'OGCounts', format='D',array = ct)
col5 = pf.Column(name = 'BGCounts', format='D',array = bkgct)

cols = pf.ColDefs([col1,col2,col3,col4,col5])
tbhdu = pf.BinTableHDU.from_columns(cols)
tbhdu.writeto(lcfile + "box2.fits",overwrite=True)
hd1.close();hd2.close();del hd1,hd2



