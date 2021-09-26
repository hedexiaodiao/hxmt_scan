#!/hxmt/soft/Develop/anaconda2/bin/python
import numpy as np
import matplotlib
matplotlib.use('Agg')
import astropy.io.fits as pf
import os,sys
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import math
#from decimal import *
#import lc_pic
#ObsID="P0101295001"
#for ObsID in range(5,6):
#	obid2="0%s"%(ObsID+1)
def reatt(attfile,Wpath,ObsID):
	OrgWpath = Wpath + '/Org'  ###
	NetWpath = Wpath + '/Net'  ###
	Net_ObsWpath = NetWpath + '/'+ObsID
	Org_attpath = "%s/%s" % (OrgWpath, 'Att')  ###
	Wcut_attfile = "%s/%s_cut_Att.fits" % (Org_attpath, ObsID)  ###perfect
	#att_all = pf.open("./%s%s/att.fits"%(ObsID,obid2))
	att_all = pf.open(attfile) 
	att_data1=att_all[1].data
	att_time1 = att_data1.field(0)
	att_ra = att_data1.field(1)
	att_dec = att_data1.field(2)
	att_dra = att_data1.field(3)
	att_ddec = att_data1.field(4)
	att_d = att_data1.field(5)
	
	att_data2=att_all[2].data
	att_time2 = att_data2.field(0)
	e_phi = att_data2.field(1)
	e_theta = att_data2.field(2)
	e_psi = att_data2.field(3)
	
	att_data3=att_all[3].data
	att_time3 = att_data3.field(0)
	q1 = att_data3.field(1)
	q2 = att_data3.field(2)
	q3 = att_data3.field(3)
	
	att_data4=att_all[4].data
	att_time4 = att_data4.field(0)
	o_x = att_data4.field(1)
	o_y = att_data4.field(2)
	o_z = att_data4.field(3)
	
	#lc_all = pf.open("./%s%s/ME/me_lc_box0_small.fits"%(ObsID,obid2))
	lc_all0 = pf.open("%s/me_lc_box0_small.fits"%(Net_ObsWpath))
	lc_all1 = pf.open("%s/me_lc_box1_small.fits"%(Net_ObsWpath))
	lc_all2 = pf.open("%s/me_lc_box2_small.fits"%(Net_ObsWpath))
	lc_data0 = lc_all0[1].data
	lc_time0 = lc_data0.field(0)
	lc_ct0 = lc_data0['Counts']
	lc_ct_e0 = lc_data0['Stat_err']
	lc_data1 = lc_all1[1].data
	lc_time1 = lc_data1.field(0)
	lc_ct1 = lc_data1.field(1)
	lc_ct_e1 = lc_data1.field(2)
	lc_data2 = lc_all2[1].data
	lc_time2 = lc_data2.field(0)
	lc_ct2 = lc_data2.field(1)
	lc_ct_e2 = lc_data2.field(2)
	mask0 = np.in1d(lc_time0,lc_time1[np.in1d(lc_time1,lc_time2)])
	mask1 = np.in1d(lc_time1,lc_time2[np.in1d(lc_time2,lc_time0)])
	mask2 = np.in1d(lc_time2,lc_time0[np.in1d(lc_time0,lc_time1)])
	#lc_time = np.array(list(set(lc_time2).union(set(lc_time1).union(set(lc_time0)))))
	#lc_time = lc_time0
	lc_time = lc_time0[mask0]
	
	col10 = pf.Column(name = 'Time', format='D',array = lc_time0[mask0])
	col20 = pf.Column(name = 'Counts', format='D',array = lc_ct0[mask0])
	col30 = pf.Column(name = 'Stat_err', format='D',array = lc_ct_e0[mask0])
	cols0 = pf.ColDefs([col10,col20,col30])
	tbhdu0 = pf.BinTableHDU.from_columns(cols0)
	lctime0 = lc_time0[mask0]
	lctime_now = lctime0[:-1]
	lctime_lag = lctime0[1:]
	lctime_mark = lctime_lag-lctime0[:-1]
	lctime_start  = lctime_lag[lctime_mark>1]
	lctime_stop  = lctime_now[lctime_mark>1]
	lcgtitstart = np.r_[lctime0[0],lctime_start]
	lcgtitstop = np.r_[lctime_stop,lctime0[-1]]
	print "gti:",lcgtitstart,lcgtitstop
	col10a = pf.Column(name = 'TSTART', format='D',array = lcgtitstart)
	col20a = pf.Column(name = 'TSTOP', format='D',array = lcgtitstop)
	cols0a = pf.ColDefs([col10a,col20a])
	prihdu = pf.PrimaryHDU()
	gtihdu = pf.BinTableHDU.from_columns(cols0a)
	tbhdulist0 = pf.HDUList([prihdu,tbhdu0,gtihdu])
	tbhdulist0.writeto("%s/me_lc_box0_small_cut.fits"%(Net_ObsWpath),clobber=True)
	
	col11 = pf.Column(name = 'Time', format='D',array = lc_time1[mask1])
	col21 = pf.Column(name = 'Counts', format='D',array = lc_ct1[mask1])
	col31 = pf.Column(name = 'Stat_err', format='D',array = lc_ct_e1[mask1])
	cols1 = pf.ColDefs([col11,col21,col31])
	tbhdu1 = pf.BinTableHDU.from_columns(cols1)
	tbhdulist1 = pf.HDUList([prihdu,tbhdu1,gtihdu])
	tbhdulist1.writeto("%s/me_lc_box1_small_cut.fits"%(Net_ObsWpath),clobber=True)
	
	col12 = pf.Column(name = 'Time', format='D',array = lc_time2[mask2])
	col22 = pf.Column(name = 'Counts', format='D',array = lc_ct2[mask2])
	col32 = pf.Column(name = 'Stat_err', format='D',array = lc_ct_e2[mask2])
	cols2 = pf.ColDefs([col12,col22,col32])
	tbhdu2 = pf.BinTableHDU.from_columns(cols2)
	tbhdulist2 = pf.HDUList([prihdu,tbhdu2,gtihdu])
	tbhdulist2.writeto("%s/me_lc_box2_small_cut.fits"%(Net_ObsWpath),clobber=True)

	
	att_time_new = lc_time
	att_ra_new = np.interp(lc_time,att_time1,att_ra)
	att_dec_new = np.interp(lc_time,att_time1,att_dec)
	att_dra_new = np.interp(lc_time,att_time1,att_dra)
	att_ddec_new = np.interp(lc_time,att_time1,att_ddec)
	att_d_new = np.interp(lc_time,att_time1,att_d)
	
	col1 = pf.Column(name = 'Time', format='D',array = lc_time)
	col2 = pf.Column(name = 'Ra', format='E',array = att_ra_new)
	col3 = pf.Column(name = 'Dec', format='E',array = att_dec_new)
	col4 = pf.Column(name = 'Delta_Ra', format='E',array = att_dra_new)
	col5 = pf.Column(name = 'Delta_Dec', format='E',array = att_ddec_new)
	col6 = pf.Column(name = 'Delta', format='E',array = att_d_new)
	cols = pf.ColDefs([col1,col2,col3,col4,col5,col6])
	tbhdu1 = pf.BinTableHDU.from_columns(cols,name='ATT_Pointing')
	
	e_phi_new = np.interp(lc_time,att_time2,e_phi)
	e_theta_new = np.interp(lc_time,att_time2,e_theta)
	e_psi_new = np.interp(lc_time,att_time2,e_psi)
	
	col12 = pf.Column(name = 'Time', format='D',array = lc_time)
	col22 = pf.Column(name = 'Euler_Phi', format='E',array = e_phi_new)
	col32 = pf.Column(name = 'Euler_Theta', format='E',array = e_theta_new)
	col42 = pf.Column(name = 'Euler_Psi', format='E',array = e_psi_new)
	cols2 = pf.ColDefs([col12,col22,col32,col42])
	tbhdu2 = pf.BinTableHDU.from_columns(cols2,name='ATT_Euler')
	
	q1_new = np.interp(lc_time,att_time3,q1)
	q2_new = np.interp(lc_time,att_time3,q2)
	q3_new = np.interp(lc_time,att_time3,q3)
	
	col13 = pf.Column(name = 'Time', format='D',array = lc_time)
	col23 = pf.Column(name = 'Q1', format='E',array = q1_new)
	col33 = pf.Column(name = 'Q2', format='E',array = q2_new)
	col43 = pf.Column(name = 'Q3', format='E',array = q3_new)
	cols3 = pf.ColDefs([col13,col23,col33,col43])
	tbhdu3 = pf.BinTableHDU.from_columns(cols3,name='ATT_Quater')
	
	o_x_new = np.interp(lc_time,att_time4,o_x)
	o_y_new = np.interp(lc_time,att_time4,o_y)
	o_z_new = np.interp(lc_time,att_time4,o_z)
	
	col14 = pf.Column(name = 'Time', format='D',array = lc_time)
	col24 = pf.Column(name = 'Omega_X', format='I',array = o_x_new)
	col34 = pf.Column(name = 'Omega_Y', format='I',array = o_y_new)
	col44 = pf.Column(name = 'Omega_Z', format='I',array = o_z_new)
	cols4 = pf.ColDefs([col14,col24,col34,col44])
	tbhdu4 = pf.BinTableHDU.from_columns(cols4,name='ATT_Omega')
	
	prihdu = pf.PrimaryHDU()
	tbhdulist = pf.HDUList([prihdu,tbhdu1,tbhdu2,tbhdu3,tbhdu4])
	#tbhdulist.writeto("./%s%s/ME/att.fits"%(ObsID,obid2),clobber=True)
	tbhdulist.writeto("%s"%(Wcut_attfile),clobber=True)


if __name__ == "__main__":
	dir_strlist = []
	with open('dir_config.txt', 'r', encoding='utf-8') as f:
		for line in f.readlines():
			line = line.strip()
			dir_strlist.append(line)
	print(dir_strlist)
	program_tree = dir_strlist[0]
	scan_tree = dir_strlist[1]
	Wpath = scan_tree  ###'/sharefs/hbkg/data/SCAN/ME'

	if len(sys.argv) < 2:
		ObsID = input("input keywords like P010129400101 :")
	elif len(sys.argv) == 2:
		ObsID = sys.argv[1]
	else:
		print('Wrong number of parameters!')
		os.exit(0)
	path = '/hxmt/work/HXMT-DATA/1L/A%s/%s/' % (ObsID[1:3], ObsID[:-5])

	###Wpath='/sharefs/hbkg/user/saina/data294/P0101294/'
	#lclist= os.listdir(Wpath)
	#lclist.sort()
	###lclist='P010129400101'
	path=Wpath+ObsID+'/ME/'
	attfile=Wpath+ObsID+"/att.fits"
	reatt(attfile,path)
