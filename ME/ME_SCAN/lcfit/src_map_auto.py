#!/hxmt/home/lick/soft/anaconda2/bin/python
'''
Copyright 2016, HXMT Ground System
All right reserved
File Name: hxmtpsf.py for sample

.......

Author: lick
Version: 0.1
Date:  2016-12-11
Email: lick@ihep.ac.cn
......
'''

import astropy.io.fits as pf
import numpy as np
from math import *
import time
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
#############read ra dec for scan stripe################


def src_map(Wpath, ObsID):
	OrgWpath = Wpath + '/Org'  ###
	NetWpath = Wpath + '/Net'  ###
	Org_obspath = "%s/%s" % (OrgWpath, ObsID)  ###
	Net_obspath = "%s/%s" % (NetWpath, ObsID)  ###
	Midd_obspath = Wpath + '/Midd/' + ObsID  ###
	Org_attpath = "%s/%s" % (OrgWpath, 'Att')  ###
	Wattfile = "%s/%s_Att.fits" % (Org_attpath, ObsID)  ###perfect
	Wcut_attfile = "%s/%s_cut_Att.fits" % (Org_attpath, ObsID)
	###srcs=np.loadtxt("%s/%s/ME/result/maxi&integral&swift.dat" % (Wpath, ObsID))
	centre_sq = np.loadtxt("%s/%s/centre.dat" % (OrgWpath, ObsID))
	num_centre = centre_sq[0]
	ra_centre = centre_sq[1]
	dec_centre = centre_sq[2]
	#print srcs,srcs.shape
	fig = plt.figure(figsize = (16,9))
	flag=0
	if os.path.exists('%s/%s_new.txt'%(Midd_obspath, ObsID)):
		src_imo = np.loadtxt("%s/%s_new.txt" % (Midd_obspath, ObsID))
		flag=1
		if src_imo.size>11:
			src_im = src_imo[-1]
		else:
			src_im = src_imo
		#print src_im
		plt.scatter(src_im[1],src_im[3],marker='*',s=50,c='crimson',edgecolor = 'crimson')
		plt.annotate('%.2f'%(src_im[9]), xy=(src_im[1],src_im[3]), xytext=(src_im[1]+0.05,src_im[3]))
	###plt.scatter(srcs[:,1],srcs[:,2],marker='o',c='seagreen',s=8,edgecolor = 'seagreen',label='srcs',alpha=0.7)
	src_cal_u = np.loadtxt("%s/%s_source_ub.txt" % (Midd_obspath, ObsID), delimiter='\t', dtype=bytes).astype(str)
	plt.scatter(src_cal_u[:,3].astype(float),src_cal_u[:,4].astype(float),marker='o',s=8,c='b',edgecolor = 'b')
	for i in range(len(src_cal_u[:,1])):
		plt.annotate('%s'%(src_cal_u[:,1][i]), xy=(src_cal_u[:,3][i].astype(float),src_cal_u[:,4][i].astype(float)), xytext=(src_cal_u[:,3][i].astype(float)+0.05,src_cal_u[:,4][i].astype(float)))
	
	if os.path.exists("%s/%s_source_bright.txt"%(Midd_obspath, ObsID)):
		src_cal_n = np.loadtxt("%s/%s_source_bright.txt" % (Midd_obspath, ObsID), delimiter='\t', dtype=bytes).astype(str)
		if src_cal_n.size<=8:
			src_cal_n=np.array([src_cal_n.tolist()])
		plt.scatter(src_cal_n[:,3].astype(float),src_cal_n[:,4].astype(float),marker='+',c='brown',s=8*(src_cal_n[:,5].astype(float)),edgecolor = 'none',label='srcs_bright')
		for i in range(len(src_cal_n[:,1])):
			plt.annotate('%s'%(src_cal_n[:,1][i]), xy=(src_cal_n[:,3][i].astype(float),src_cal_n[:,4][i].astype(float)), xytext=(src_cal_n[:,3][i].astype(float)+0.05,src_cal_n[:,4][i].astype(float)))
	#for j in range(len(src_im[:,0])):
		#plt.annotate('%.2f'%src_im[:,4][j], xy=(src_im[:,0][j],src_im[:,2][j]), xytext=(src_im[:,0][j]+0.2,src_im[:,2][j]-0.5))
	
	att = pf.open(Wcut_attfile)
	att_o = pf.open(Wattfile)
	att_ra=att[1].data.field(1)
	att_dec=att[1].data.field(2)
	att_ora=att_o[1].data.field(1)
	att_odec=att_o[1].data.field(2)
	ramin = 0
	decmin = 0
	part_ra=[]
	part_dec=[]
	k=0
	plt.plot(att_ora,att_odec,'y:',alpha=0.5)
	plt.scatter(att_ora[0],att_odec[0],c='y',s=15,marker='>',alpha=0.7)
	for i in range(len(att_ra)-1):
		part_ra.append(att_ra[i])
		part_dec.append(att_dec[i])
		if sqrt((att_ra[i]-att_ra[i+1])**2+(att_dec[i]-att_dec[i+1])**2)>=0.5:
			plt.plot(part_ra,part_dec,'r-')
			k=0
			part_ra=[]
			part_dec=[]
			continue
		if i+1 == len(att_ra)-1:
			plt.plot(part_ra,part_dec,'r-')
	
	print ObsID, ra_centre, dec_centre
	plt.title(ObsID, loc="left")
	plt.xlabel("R.A., degree")
	plt.ylabel("Dec., degree")
	#plt.legend(bbox_to_anchor=(1.12, 1.0))
	plt.grid()
	#plt.show()
	plt.savefig('%s/%s_smap.png' % (Midd_obspath, ObsID), dpi = 400)
	plt.close('all')


if __name__ == "__main__":
	Wpath='/sharefs/hbkg/user/saina/data294/P0211007/'
	lclist= os.listdir(Wpath)
	lclist.sort()
	for i in lclist:
		obid=i
		if os.path.exists('%s/%s/ME/result/%s_smap.png'%(Wpath,obid,obid)):
			continue
		if os.path.exists('%s/%s/ME/result/%s_ub_chi2.dat'%(Wpath,obid,obid)):
			src_map(Wpath,obid)
