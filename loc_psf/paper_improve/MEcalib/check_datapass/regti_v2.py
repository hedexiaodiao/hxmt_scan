#!/hxmt/soft/Develop/anaconda2/bin/python
import numpy as np
import astropy.io.fits as pf
import math
import glob
import delect_motionless
import calc_intersection

#ObsID = "P0101295002"
#Wpath = "/hxmt/work/HXMT-DATA/1R/PIPROD/A01/P0101295/"+ObsID+"/"
#for i in range(8):
def regti(meevtfile, Wpath, ObsID):
	OrgWpath = Wpath + '/Org'###
	Org_attpath = "%s/%s" % (OrgWpath, 'Att')
	Org_gtipath = "%s/%s" % (OrgWpath, 'GTI')  ###
	Wgtifile = "%s/%s_gtiv2.fits" % (Org_gtipath, ObsID)  ###
	Wme_gtifile = "%s/%s_me_gti.fits" % (Org_gtipath, ObsID)  ###perfect
	Wattfile = "%s/%s_Att.fits" % (Org_attpath, ObsID)

	#obid2 = "0%s"%(i+1)
	#meevtfile_list = glob.glob(r'%s/%s%s*/*/*ME-Evt*V1*'%(Wpath,ObsID,obid2))
	#for i in meevtfile_list:
	#	meevtfile = i
	#print meevtfile
	evt_all = pf.open(meevtfile)
	soft_gti = pf.open(Wgtifile)
	lenfile = len(soft_gti)
	softhdu = soft_gti[lenfile-1]
	evt_data = evt_all[1].data
	evt_time = evt_data.field(0)
	start_t = min(evt_time)
	start_t=np.array([start_t])
	stop_t = max(evt_time)
	stop_t=np.array([stop_t])
	print start_t,stop_t
	motion_start, motion_stop = delect_motionless.motionless(Wattfile)

	new_start_t,new_stop_t = calc_intersection.calc_intersection(start_t,stop_t,motion_start,motion_stop)

	tbhdu=[[] for j in range(lenfile-2)]
	for i in range(1,lenfile-1):
		col1 = pf.Column(name = 'START', format='D',array = new_start_t)
		col2 = pf.Column(name = 'STOP', format='D',array = new_stop_t)
		cols = pf.ColDefs([col1,col2])
		tbhdu[i-1] = pf.BinTableHDU.from_columns(cols,name = 'GTI%s'%(i-1))
	
	prihdu = pf.PrimaryHDU()
	list1=[prihdu]
	for k in range(len(tbhdu)):
		list1=np.append(list1,tbhdu[k])
	
	list1=np.append(list1,softhdu)
	list1=list(list1)
	tbhdulist = pf.HDUList(list1)
	#tbhdulist.writeto("./%s%s/ME/me_gti.fits"%(ObsID,obid2),clobber=True)
	tbhdulist.writeto(Wme_gtifile,clobber=True)

