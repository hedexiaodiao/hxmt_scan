#!/hxmt/soft/Develop/anaconda2/bin/python
import numpy as np
import astropy.io.fits as pf
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.collections as collections
import math
import astropy.io.fits as pf
import sys
import xlrd
import os

def gen_lc_me(ehkfile, gtifile, Wpath, ObsID):
	Program_dir = '/sharefs/hbkg/user/luoqi/HXMT_SCAN/ME/ME_SCAN'
	OrgWpath = Wpath + '/Org'  ###
	NetWpath = Wpath + '/Net'  ###
	Net_ObsWpath = NetWpath + '/' + ObsID
	Org_obspath = "%s/%s" % (OrgWpath, ObsID)  ###
	Org_screenpath = "%s/%s" % (OrgWpath, 'Screen')  ###
	Wbaddetfile = "%s/%s_me_bdet.fits" % (Org_obspath, ObsID)  ###bad det

	timelim_fits = pf.open(gtifile)
	timelim = timelim_fits[1].data
	tstart = timelim.field(0)
	tstop = timelim.field(1)
	ehk_fits = pf.open(ehkfile)
	ehk_tb = ehk_fits[1].data
	ehk_tb = ehk_tb[(ehk_tb.field(0)>=tstart)&(ehk_tb.field(0)<tstop)]
	ehk_time = ehk_tb.field(0)
	#ehk_time=ehk_time-min(ehk_time)

	###&(ehk_tb.field("SUN_ANG")>=10)
	###&(ehk_tb.field("SAA_FLAG")==0)
	###&(ehk_tb.field("MOON_ANG")>=4)
	mask = (ehk_tb.field("ELV")>5)&(ehk_tb.field("T_SAA")>=200)&(ehk_tb.field("TN_SAA")>=100)\
	&(ehk_tb.field("COR2")>=8)\
	&(ehk_tb.field("ANG_DIST")<=359)&\
	((ehk_tb.field("SAT_LAT")<31)|(ehk_tb.field("SAT_LAT")>38))\
	&((ehk_tb.field("SAT_LON")>245)|(ehk_tb.field("SAT_LON")<288))\
	&(ehk_tb.field("SAT_LAT")>=-36.5)&(ehk_tb.field("SAT_LAT")<=36.5)
	#&(ehk_tb.field("DYE_ELV")>=0)
	gti_all = ehk_time[mask]
	gti_now = gti_all[:-1]
	gti_lag = gti_all[1:]
	gti_mark = gti_lag-gti_all[:-1]
	gti_start  = gti_lag[gti_mark>1]
	gti_stop  = gti_now[gti_mark>1]
	#print gti_all.shape,gti_lag.shape,gti_mark.shape,gti_mark,gti_start,gti_stop
	gtitstart = np.r_[gti_all[0],gti_start]
	gtitstop = np.r_[gti_stop,gti_all[-1]]
	#print gtitstart,gtitstop
	fig = plt.figure(figsize = (20,10))
	gs = gridspec.GridSpec(1, 1,bottom=0.05,top=0.96,left=0.05,right=0.95,hspace=0.05)
	gs0 = gridspec.GridSpecFromSubplotSpec(6, 1,subplot_spec=gs[0], hspace=0.)
	ax00 = plt.subplot(gs0[0])
	ax01 = plt.subplot(gs0[1],sharex=ax00)
	ax02 = plt.subplot(gs0[2],sharex=ax00)
	ax03 = plt.subplot(gs0[3],sharex=ax00)
	ax04 = plt.subplot(gs0[4],sharex=ax00)
	ax05 = plt.subplot(gs0[5],sharex=ax00)
	ax01.plot(ehk_time,ehk_tb.field("ELV"),label='ELV')
	ax01.plot(ehk_time[ehk_tb.field("ELV")>5],ehk_tb.field("ELV")[ehk_tb.field("ELV")>5],'.')
	ax01.legend()
	ax02.plot(ehk_time,ehk_tb.field("SAT_LAT"),label='SAT_LAT')
	ax02.plot(ehk_time[((ehk_tb.field("SAT_LAT")<31)|(ehk_tb.field("SAT_LAT")>38))&((ehk_tb.field("SAT_LAT")>=-36.5)&(ehk_tb.field("SAT_LAT")<=36.5))],ehk_tb.field("SAT_LAT")[((ehk_tb.field("SAT_LAT")<31)|(ehk_tb.field("SAT_LAT")>38))&((ehk_tb.field("SAT_LAT")>=-36.5)&(ehk_tb.field("SAT_LAT")<=36.5))],'.')
	ax02.legend()
	ax03.plot(ehk_time,ehk_tb.field("SAT_LON"),label='SAT_LON')
	ax03.plot(ehk_time[((ehk_tb.field("SAT_LON")>245)|(ehk_tb.field("SAT_LON")<228))],ehk_tb.field("SAT_LON")[((ehk_tb.field("SAT_LON")>245)|(ehk_tb.field("SAT_LON")<228))],'.')
	ax03.legend()
	ax04.plot(ehk_time,ehk_tb.field("COR2"),label='COR')
	ax04.plot(ehk_time[ehk_tb.field("COR2")>=8],ehk_tb.field("COR2")[ehk_tb.field("COR2")>=8],'.')
	ax04.legend()
	ax05.plot(ehk_time,ehk_tb.field("SAA_FLAG"),label='SAA_FLAG')
	ax05.legend()
	
	###################################################################################
	if os.path.exists('%s/%s_hptime_rebin.txt'%(Org_screenpath, ObsID)):###
		bdgti=np.loadtxt('%s/%s_hptime_rebin.txt' % (Org_screenpath, ObsID))###
		bdstr=bdgti[0]
		bdstp=bdgti[1]
		if bdstr.size>1:
			numbdg=len(bdgti[0])
		else:
			bdstr=np.array([bdgti[0]])
			bdstp=np.array([bdgti[1]])
			numbdg=1
		gtitstart0=gtitstart
		gtistop0=gtitstop
		print bdgti,gtitstart
		for i in range(numbdg):
			print 'turn ',i
			numgti=gtitstart.size
			for j in range(numgti):
				if (gtitstart[j]>=bdstp[i])|(gtitstop[j]<=bdstr[i]):
					#print j,' j1--exclude'
					continue
				if gtitstart[j]>=bdstr[i]:
					gtitstart[j]=bdstp[i]
					print j,' j2--last'
					continue
				if gtitstop[j]<=bdstp[i]:
					gtitstop[j]=bdstr[i]
					print j,' j3--first'
				else:
					gtitstop=np.append(gtitstop,bdstr[i])
					gtitstart=np.append(gtitstart,bdstp[i])
					gtitstop.sort()
					gtitstart.sort()
					print j,' j4--middle and add'
	
	############################################################################################
	data = xlrd.open_workbook('%s/Blank_sky/pixel_class_mu.xlsx'%(Program_dir))###
	classall=data.sheets()[3]
	numpix=classall.col_values(0)
	classnum=classall.col_values(3)
	ptype = np.loadtxt('%s/Blank_sky/pixel_type.txt'%(Program_dir))
	uup=ptype[ptype[:,1]>1][:,0]
	class7=np.array(numpix)[(np.array(classnum)>6)]
	allnum=np.array(range(1728))
	mask1=np.in1d(allnum,uup)
	mask2=np.in1d(allnum,class7)
	sfile=pf.open(Wbaddetfile)###bad det file
	bdnum=sfile[1].data.field(0)
	nup=allnum[mask1|mask2]
	#mask3=np.in1d(bdnum,nup)
	mask4=np.in1d(allnum,bdnum,invert=True)
	usp=allnum[np.in1d(allnum,nup,invert=True)]
	gp=allnum[mask4]
	gp2=gp[np.in1d(gp,nup,invert=True)]
	r0=float(len(usp[usp<576]))/float(len(gp2[gp2<576]))
	r1=float(len(usp[(usp>=576)&(usp<1152)]))/float(len(gp2[(gp2>=576)&(gp2<1152)]))
	r2=float(len(usp[usp>=1152]))/float(len(gp2[gp2>=1152]))
	rall=[r0,r1,r2]
	for i in range(3):
		lcfits = pf.open("%s/me_lc_box%s_small_bkg.fits" % (Net_ObsWpath, i))
		lctime_grade_small_short = lcfits[1].data.field(0)
		#print "lc:",i,lctime_grade_small_short.shape,lctime_grade_small_short,lctime_grade_small_short[0],lctime_grade_small_short[-1]
		lc_grade_small_short = lcfits[1].data.field(1)
		lc_grade_small_short_e = lcfits[1].data.field(2)
		lc_grade_small = []
		lctime_grade_small = []
		lc_grade_small_e = []
		for k in range(len(gtitstart)):
			lctime_grade_small = np.r_[lctime_grade_small,lctime_grade_small_short[(lctime_grade_small_short>=gtitstart[k])&(lctime_grade_small_short<=gtitstop[k])]]
			lc_grade_small = np.r_[lc_grade_small,lc_grade_small_short[(lctime_grade_small_short>=gtitstart[k])&(lctime_grade_small_short<=gtitstop[k])]]
			lc_grade_small_e = np.r_[lc_grade_small_e,lc_grade_small_short_e[(lctime_grade_small_short>=gtitstart[k])&(lctime_grade_small_short<=gtitstop[k])]]
		lctime_grade_small0 = lctime_grade_small-min(ehk_time)
		lctime_grade_small_short0=lctime_grade_small_short-min(ehk_time)
		#print "new:",lctime_grade_small.shape,lctime_grade_small
		ax00.plot(lctime_grade_small_short,lc_grade_small_short*rall[i],label='Box%s'%i)
		
		col1 = pf.Column(name = 'Time', format='D',array = lctime_grade_small)
		col2 = pf.Column(name = 'Counts', format='D',array = lc_grade_small*rall[i])
		col3 = pf.Column(name = 'Stat_err', format='D',array = lc_grade_small_e*rall[i])
		cols = pf.ColDefs([col1,col2,col3])
		tbhdu = pf.BinTableHDU.from_columns(cols)
		###########gti##############################################
		lctime = lctime_grade_small
		lctime_now = lctime[:-1]
		lc = lc_grade_small
		lctime_lag = lctime[1:]
		lctime_mark = lctime_lag-lctime[:-1]
		lctime_start  = lctime_lag[lctime_mark>1]
		lctime_stop  = lctime_now[lctime_mark>1]
		#print lctime.shape,lctime_lag.shape,lctime_mark.shape,lctime_mark,lctime_start,lctime_stop
		lcgtitstart = np.r_[lctime[0],lctime_start]
		lcgtitstop = np.r_[lctime_stop,lctime[-1]]
		print "gti:",lcgtitstart,lcgtitstop
		col12 = pf.Column(name = 'TSTART', format='D',array = lcgtitstart)
		col22 = pf.Column(name = 'TSTOP', format='D',array = lcgtitstop)
		cols2 = pf.ColDefs([col12,col22])
		prihdu = pf.PrimaryHDU()
		gtihdu = pf.BinTableHDU.from_columns(cols2)
		tbhdulist = pf.HDUList([prihdu,tbhdu,gtihdu])
		tbhdulist.writeto("%s/me_lc_box%s_small.fits" % (Net_ObsWpath, i), clobber=True)
	
	#ax00.plot(ehk_time,np.array(mask,dtype=int)*(-15),c='k',lw=0.7,label='GTI')
	#ax00.legend()
	head_file=pf.open("%s/me_0_small_g0_0-17.lc" % (Net_ObsWpath))
	head_or=head_file[1].header
	date=head_or['DATE-OBS']
	collection = collections.BrokenBarHCollection.span_where(ehk_time, ymin=-100, ymax=3000, where=np.array(mask,dtype=int) > 0, alpha=0.5)
	ax00.add_collection(collection)
	if os.path.exists('%s/%s_hptime_rebin.txt'%(Org_screenpath, ObsID)):###
		mark1=(ehk_time>=bdstr[0])&(ehk_time<=bdstp[0])
		for k in range(1,len(bdstr)):
			mark1=mark1|(ehk_time>=bdstr[k])&(ehk_time<=bdstp[k])
		collection1 = collections.BrokenBarHCollection.span_where(ehk_time, ymin=-100, ymax=3000, where=np.array(mark1,dtype=int) > 0, alpha=0.3,color='orange')
		ax00.add_collection(collection1)
	ax00.legend()
	ax00.set_ylim(min(lc_grade_small_short*rall[0])-20,max(lc_grade_small_short*rall[0])+30)
	ax00.set_title('%s    %s' % (ObsID, date))
	plt.xlim(min(ehk_time),max(ehk_time))
	ax05.set_xlabel('Time(s)',fontsize=12)
	ax00.tick_params(labelbottom=False)
	ax01.tick_params(labelbottom=False)
	ax02.tick_params(labelbottom=False)
	ax03.tick_params(labelbottom=False)
	ax04.tick_params(labelbottom=False)
	ax00.grid()
	ax01.grid()
	ax02.grid()
	ax03.grid()
	ax04.grid()
	ax05.grid()
	plt.savefig('%s/me_small_all_%s.png' % (Net_ObsWpath, ObsID), dpi=400)
	plt.close('all')

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
	#path = '/hxmt/work/HXMT-DATA/1L/A%s/%s/' % (ObsID[1:3], ObsID[:-5])

	###Wpath='/sharefs/hbkg/user/saina/data294/P0101294/'
	#lclist= os.listdir(Wpath)
	#lclist.sort()
	###lclist='P010129400101'
	path=Wpath+ObsID+'/ME/'###
	ehkfile=Wpath+ObsID+"/ehk.fits"
	gtifile=path+'me_gti.fits'
	gen_lc_me(ehkfile,gtifile,path,ObsID)