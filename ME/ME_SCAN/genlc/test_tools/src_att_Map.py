import sys
sys.path.append("/hxmt/home/saina/tar/test_tools/")
import astropy.io.fits as pf
import numpy as np
from math import *
import time
import Quat as quat
import matplotlib.pyplot as plt
#############read ra dec for scan stripe################
def src_map(infile,ft,slist):
	hdulist = pf.open(infile)
	tb1 = hdulist[1].data
	tb3 = hdulist[3].data
	qtime = tb3.field(0)
	q1_list = tb3.field(1)
	q2_list = tb3.field(2)
	q3_list = tb3.field(3)
	stripe_ra_list = []
	stripe_dec_list = []
	print "q list size: ",q1_list.shape
	fatt=open("att.dat",'w')
	fq=open("q.dat",'w')
	for i in range(0,q1_list.shape[0]):
		qu = [q1_list[i],q2_list[i],q3_list[i]]
		qu = qu+[sqrt(fabs(1-qu[0]*qu[0]-qu[1]*qu[1]-qu[2]*qu[2]))]
		quat1 = quat.Quat(quat.normalize(qu))
		stripe_ra_list.append(quat1.ra)
		stripe_dec_list.append(quat1.dec)
	time.sleep(5)
	fatt.close()
	fq.close()
	stripe_ra_list = np.array(stripe_ra_list)
	stripe_dec_list = np.array(stripe_dec_list)
	ra_centre = stripe_ra_list.mean()
	dec_centre = stripe_dec_list.mean()
	ralb = stripe_ra_list.min()
	rahb = stripe_ra_list.max()
	declb = stripe_dec_list.min()
	dechb = stripe_dec_list.max()
	#src_list = pf.open("/hxmt/home/saina/tar/test_tools/Srcs_IGR_SWIFT_201811.fits")
	src_list = pf.open(slist)
	src_map = src_list[1].data
	src_name = src_map.field(0)
	src_ra = src_map.field(1)
	src_dec = src_map.field(2)
	src_flux = src_map.field(6)
	Yary=[]
	num_list=[]
	for i,itra in enumerate(src_ra):
		if ((((itra - stripe_ra_list)**2 + (src_dec[i]-stripe_dec_list)**2)<=17)).sum()>0:
			Yary.append(1)
			num_list.append(i)
		else:
			Yary.append(0)
	print sum(Yary)
	Yary=(np.array(Yary)>0)
	#src_num = np.nonzero(Yary)[0]
	src_ra_sel = src_ra[Yary&(src_flux>=ft)]
	src_dec_sel = src_dec[Yary&(src_flux>=ft)]
	src_flux_sel = src_flux[Yary&(src_flux>=ft)]
	src_name_sel = src_name[Yary&(src_flux>=ft)]
	src_flux_tot = src_flux_sel
	src_map_sel = np.vstack((src_flux_tot,src_ra_sel,src_dec_sel))
	return src_map_sel.T,stripe_ra_list,stripe_dec_list,src_name_sel

if __name__ == "__main__":
	obid='P0101294126'
	path="/hxmt/work/USERS/saina/P0101294/"+obid+"/"
	flux_thres = 5.0
	srcs,ralist,declist,src_name = src_map('%s/ME/att.fits'%path,flux_thres)
	print "srcs:::::::",srcs
	plt.scatter(srcs[:,1],srcs[:,2],marker='o',s=10,c = 'y',edgecolor = 'k')
	att = pf.open('%s/ME/att.fits'%path)
	ra=att[1].data.field(1)
	dec=att[1].data.field(2)
	sel=[True]*len(ra)
	for i in srcs:
		sel_o=(((i[1]-ra)**2+(i[2]-dec)**2)>17)
		sel=sel&sel_o
	
	plt.plot(ra,dec,'r-')
	plt.scatter(ra[sel],dec[sel],s=1,c='b')
	#plt.scatter(ralist,declist,c='b',s=3)
	plt.grid()
	plt.show()
