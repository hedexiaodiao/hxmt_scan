#!/usr/bin/env python
#coding:utf-8
import numpy as np
import xlrd
import astropy.io.fits as pf
###-----------------------------------------------
'''
归一化一开始对不上是因为没有剔除放射源
'''


ObsID = 'P020204131503'
Program_dir = '/sharefs/hbkg/user/luoqi/HXMT_SCAN/ME/ME_SCAN'
Org_obspath = '/sharefs/hbkg/data/SCAN/luoqi/calib/7_20/Org/%s'%(ObsID)

Wbaddetfile = "%s/%s_me_bdet.fits" % (Org_obspath, ObsID)

data = xlrd.open_workbook('%s/Blank_sky/pixel_class_mu.xlsx'%(Program_dir))###
classall=data.sheets()[3]
numpix=classall.col_values(0)
classnum=classall.col_values(3)
ptype = np.loadtxt('%s/Blank_sky/pixel_type.txt'%(Program_dir))

print('ptype:',ptype)
uup=ptype[ptype[:,1]>1][:,0]
print('uup:',uup)

class7=np.array(numpix)[(np.array(classnum)>6)]
allnum=np.array(range(1728))
mask1=np.in1d(allnum,uup)
mask2=np.in1d(allnum,class7)
print(allnum[mask1],allnum[mask2],len(allnum[mask1]),len(allnum[mask2]))

sfile=pf.open(Wbaddetfile)###bad det file
bdnum=sfile[1].data.field(0)
nup=allnum[mask1|mask2]###符合两个表格里条件之一的，都是坏探测器
print('nup:',nup,len(nup))

#mask3=np.in1d(bdnum,nup)
mask4=np.in1d(allnum,bdnum,invert=True)###坏探测器不在所有探测数里
usp=allnum[np.in1d(allnum,nup,invert=True)]###除去两个表格条件的好探测器
###print('mask4',allnum[mask4])

gp=allnum[mask4]###所有去除坏探测器fits筛选的好探测器
gp2=gp[np.in1d(gp,nup,invert=True)]###fits筛选后再经历两个表格筛选，剩余的好探测器
print('gp, is mask4:',gp)

print('gp2:',gp2,len(gp2))
print('gp2>=576 <152 box1,', len(gp2[(gp2>=576)&(gp2<1152)]))

r0=float(len(usp[usp<576]))/float(len(gp2[gp2<576]))
r1=float(len(usp[(usp>=576)&(usp<1152)]))/float(len(gp2[(gp2>=576)&(gp2<1152)]))
r2=float(len(usp[usp>=1152]))/float(len(gp2[gp2>=1152]))
rall=[r0,r1,r2]

###print('ptype:',ptype)
print('uup:',len(uup[(uup>=576)&(uup<1152)]))
dmask1 = allnum[mask1]
dmask2 = allnum[mask2]
box1usp = usp[(usp>=576)&(usp<1152)]
box1gp2 = gp2[(gp2>=576)&(gp2<1152)]

print('mask1 mask2 num,:',len(dmask1[(dmask1>=576)&(dmask1<1152)]),len(dmask2[(dmask2>=576)&(dmask2<1152)]))
print('uup box1,',uup[(uup>=576)&(uup<1152)])
print('nup num:',len(nup[(nup>=576)&(nup<1152)&(nup<=831)&(nup>=928)]))
print('gp num, is mask4:',len(gp[(gp>=576)&(gp<1152)]))
print('gp2>=576 <152 box1,', len((gp2[((gp2>=576)&(gp2<=831))|((gp2<1152)&(gp2>=928))])))
print('usp,',len(usp[((usp>=576)&(usp<=831))|((usp<1152)&(usp>=928))]))