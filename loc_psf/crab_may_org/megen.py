import sys
import numpy as np
import xlrd
import os
import glob
from astropy.io import fits as pf
#bdfile = '/home/hxmt/hxmtsoft2/hxmtsoftv2.02/install/x86_64-pc-linux-gnu-libc2.12/refdata/medetectorstatus.fits'
from tool import *

obsid = sys.argv[1]
pfile = "/home/hxmt/nangyi/pfiles/me%s"%str(obsid)
os.system('mkdir %s'%pfile);os.system("rm %s/*"%pfile);
pfilepath = "%s;/home/hxmt/hxmtsoft2/hxmtsoftv2.02/install/x86_64-pc-linux-gnu-libc2.12/syspfiles"%pfile
cmdinit = 'sleep 1;source /sharefs/hbkg/user/nangyi/hxmtsoft_v2.02.sh ;export PFILES="%s" '%(pfilepath)

if obsid[:11]=='P0101299011':bdfile='/sharefs/hbkg/user/saina/data295/P0101295/P010129503801/bd_P010129503801.fits'
if obsid[:11]=='P0101299012':bdfile='/sharefs/hbkg/user/saina/data295/P0101295/P010129503901/bd_P010129503901.fits'
if obsid[:11]=='P0101299013':bdfile='/sharefs/hbkg/user/saina/data295/P0101295/P010129506501/bd_P010129506501.fits'

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

cmd=''
if not os.path.exists(pifile):
    cmd='; mepical evtfile=%s tempfile=%s outfile=%s'%(basefile['EVT'],basefile['TH'],pifile)

if not os.path.exists(grdfile):
    cmd +=' ;megrade evtfile=%s outfile=%s deadfile=%s binsize=1'%(pifile,grdfile,deadfile)

if not os.path.exists(gtifile):
    cmd +=' ;megtigen tempfile="%s" ehkfile="%s" outfile="%s" defaultexpr=NONE expr="ELV>10&&COR>8&&SAA_FLAG==0&&T_SAA>300&&TN_SAA>300&&ANG_DIST<=0.05" clobber=yes history=yes'%(basefile['TH'], basefile['EHK'], gtifile)
    cmd +=' ;megti %s %s %s $HEADAS/refdata/medetectorstatus.fits %s'%(grdfile, gtifile, newgtifile, newbdfile)

#newbdfile = bdfile
if not os.path.exists(screenfile):
    cmd +=' ;mescreen evtfile=%s gtifile=%s outfile=%s baddetfile=%s userdetid="0-53"'%(grdfile,newgtifile,screenfile,newbdfile)

os.system(cmdinit + cmd)
thd = pf.open(screenfile)#,endcard=False)
try:
    print len(thd)
except:
    pass

startm = int(thd[2].data.field(0)[0])+1
stoptm = int(thd[2].data.field(1)[-1])-1

thd.close();del thd

cmd = '; melcgen evtfile=%s deadfile=%s outfile=%s userdetid="0-7,11-17;18-25,29-35;36-43,47-53" starttime=%s stoptime=%s minPI=68 maxPI=631 binsize=1 deadcorr=yes'%(screenfile,deadfile,lcfile,startm,stoptm)
cmd += ' ; ls -t %s_g0*.lc | head -1 >  %s'%(lcfile,lcnametxt)
cmd += ' ; python mebkgmap_2012005_lc_err.py lc %s %s %s %s %s %s 68 631 %s %s'%(screenfile,basefile['EHK'],newgtifile,deadfile,basefile['TH'],lcnametxt,bkgfile,newbdfile)
os.system(cmdinit + cmd)
allnum=np.array(range(1728))
data = xlrd.open_workbook('pixel_class_mu.xlsx')
classall=data.sheets()[3]
numpix=classall.col_values(0)
classnum=classall.col_values(3)
class7=np.array(numpix)[(np.array(classnum)>6)]
ptype = np.loadtxt('pixel_type.txt')
uup=ptype[ptype[:,1]>1][:,0]
mask1=np.in1d(allnum,uup)
mask2=np.in1d(allnum,class7)
sfile=pf.open("%s"%newbdfile)
bdnum=sfile[1].data.field(0)
nup=allnum[mask1|mask2]
mask4=np.in1d(allnum,bdnum,invert=True)
usp=allnum[np.in1d(allnum,nup,invert=True)]
gp=allnum[mask4]
gp2=gp[np.in1d(gp,nup,invert=True)]
r0=float(len(usp[usp<576]))/float(len(gp2[gp2<576]))
r1=float(len(usp[(usp>=576)&(usp<1152)]))/float(len(gp2[(gp2>=576)&(gp2<1152)]))
r2=float(len(usp[usp>=1152]))/float(len(gp2[gp2>=1152]))
rall=[r0,r1,r2]
print rall
#r0,r1,r2=1,1,1
hd1=pf.open('%s_g0_0-17.lc'%lcfile)
hd2=pf.open('%s_box0.lc'%bkgfile)
yp = (hd1[1].data['FRACEXP']==1.)
tm = hd1[1].data['TIME'][yp]
ct = hd1[1].data['COUNTS'][yp]
cter = hd1[1].data['ERROR'][yp]
bkgct = hd2[1].data['RATE'][yp]
bkger = hd2[1].data['Error'][yp]
netct = ct-bkgct
neter = np.sqrt(cter**2+bkger**2)
col1 = pf.Column(name = 'Time', format='D',array = tm)
col2 = pf.Column(name = 'Counts', format='D',array = netct*r0)
col3 = pf.Column(name = 'Stat_err', format='D',array = neter*r0)
col4 = pf.Column(name = 'Cts_org', format='D',array = ct)
col5 = pf.Column(name = 'Cts_err', format='D',array = cter)
col6 = pf.Column(name = 'Bkg', format='D',array = bkgct)
col7 = pf.Column(name = 'Bkg_err', format='D',array = bkger)
cols = pf.ColDefs([col1,col2,col3,col4,col5,col6,col7])
tbhdu = pf.BinTableHDU.from_columns(cols)
tbhdu.writeto(lcfile + "box0.fits",overwrite=True)
hd1.close();hd2.close();del hd1,hd2


hd1=pf.open('%s_g1_18-35.lc'%lcfile)
hd2=pf.open('%s_box1.lc'%bkgfile)
tm = hd1[1].data['TIME'][yp]
ct = hd1[1].data['COUNTS'][yp]
cter = hd1[1].data['ERROR'][yp]
bkgct = hd2[1].data['RATE'][yp]
bkger = hd2[1].data['Error'][yp]
netct = ct-bkgct
neter = np.sqrt(cter**2+bkger**2)
col1 = pf.Column(name = 'Time', format='D',array = tm)
col2 = pf.Column(name = 'Counts', format='D',array = netct*r1)
col3 = pf.Column(name = 'Stat_err', format='D',array = neter*r1)
col4 = pf.Column(name = 'Cts_org', format='D',array = ct)
col5 = pf.Column(name = 'Cts_err', format='D',array = cter)
col6 = pf.Column(name = 'Bkg', format='D',array = bkgct)
col7 = pf.Column(name = 'Bkg_err', format='D',array = bkger)
cols = pf.ColDefs([col1,col2,col3,col4,col5,col6,col7])
tbhdu = pf.BinTableHDU.from_columns(cols)
tbhdu.writeto(lcfile + "box1.fits",overwrite=True)
hd1.close();hd2.close();del hd1,hd2


hd1=pf.open('%s_g2_36-53.lc'%lcfile)
hd2=pf.open('%s_box2.lc'%bkgfile)
tm = hd1[1].data['TIME'][yp]
ct = hd1[1].data['COUNTS'][yp]
cter = hd1[1].data['ERROR'][yp]
bkgct = hd2[1].data['RATE'][yp]
bkger = hd2[1].data['Error'][yp]
netct = ct-bkgct
neter = np.sqrt(cter**2+bkger**2)
col1 = pf.Column(name = 'Time', format='D',array = tm)
col2 = pf.Column(name = 'Counts', format='D',array = netct*r2)
col3 = pf.Column(name = 'Stat_err', format='D',array = neter*r2)
col4 = pf.Column(name = 'Cts_org', format='D',array = ct)
col5 = pf.Column(name = 'Cts_err', format='D',array = cter)
col6 = pf.Column(name = 'Bkg', format='D',array = bkgct)
col7 = pf.Column(name = 'Bkg_err', format='D',array = bkger)
cols = pf.ColDefs([col1,col2,col3,col4,col5,col6,col7])
tbhdu = pf.BinTableHDU.from_columns(cols)
tbhdu.writeto(lcfile + "box2.fits",overwrite=True)
hd1.close();hd2.close();del hd1,hd2



