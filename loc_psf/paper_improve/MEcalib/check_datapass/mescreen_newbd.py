import numpy as np
import astropy.io.fits as pf
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import xlrd
import os
from math import *
import time
import sys
import matplotlib.gridspec as gridspec
import matplotlib.collections as collections
sys.path.append("/sharefs/hbkg/user/luoqi/HXMT_SCAN/ME/ME_SCAN/genlc/test_tools/")
###sys.Wpath.append("/hxmt/home/saina/tar/test_tools/")
import Quat as quat
import src_att_Map
#######################################class & function(excel)####################################


def bkg_gti(Wpath, ObsID, ft,Program_dir):
	OrgWpath = Wpath + '/Org'  ###
	Org_attpath = "%s/%s" % (OrgWpath, 'Att')  ###
	Org_gtipath = "%s/%s" % (OrgWpath, 'GTI')  ###
	Wattfile = "%s/%s_Att.fits" % (Org_attpath, ObsID)  ###perfect
	Wgtifile = "%s/%s_gtiv2.fits" % (Org_gtipath, ObsID)  ###
	Wbkg_gtifile = "%s/%s_me_gti_bkg.fits" % (Org_gtipath, ObsID)  ###

	attfile=pf.open(Wattfile)
	att=attfile[1].data
	atttime=att.field(0)
	#sra=83.633
	#sdec=22.015
	attra=att.field(1)
	attdec=att.field(2)
	#slist='/sharefs/hbkg/user/saina/tar/psf_201904/Srcs_IGR_SWIFT_06.fits'
	slist='%s/LE_SOURCE_LIST.FITS'%(Program_dir)
	path= Wpath +'/' + ObsID###+"/"
	srcs,ralist,declist,src_name = src_att_Map.src_map(Wattfile,ft,slist)
	print "srcs:::::::",srcs
	sel=np.array([True]*len(attra))
	if len(srcs)!=0:
		for i in srcs:
			sel_o=(((i[1]-attra)**2+(i[2]-attdec)**2)>17)
			sel=sel&sel_o
	
	atttime_now = atttime[sel][:-1]
	atttime_lag = atttime[sel][1:]
	atttime_mark = atttime_lag-atttime[sel][:-1]
	if (atttime[sel].size>0):
		atttime_start = np.r_[atttime[sel][0],atttime_lag[atttime_mark>1]]
		atttime_stop  = np.r_[atttime_now[atttime_mark>1],atttime[sel][-1]]
	else:
		atttime_start=atttime[0]
		atttime_stop=atttime[0]
		print ObsID, ' no att time'
	sgti=pf.open(Wgtifile)
	lenfile = len(sgti)
	softhdu = sgti[lenfile-1]
	sgti_str=sgti[1].data.field(0)
	sgti_stp=sgti[1].data.field(1)
	gtistr=np.r_[sgti_str,atttime_start]
	gtistp=np.r_[sgti_stp,atttime_stop]
	new_str=np.array([])
	new_stp=np.array([])
	for j in range(len(sgti_str)):
		for i in range(len(atttime_start)):
			if (sgti_stp[j]<atttime_start[i])|(sgti_str[j]>atttime_stop[i]):
				continue
			if (atttime_start[i]>=sgti_str[j])&(atttime_start[i]<=sgti_stp[j]):
				if (atttime_stop[i]<sgti_stp[j]):
					new_str=np.append(new_str,atttime_start[i])
					new_stp=np.append(new_stp,atttime_stop[i])
					print 1,": att in"
				else:
					new_str=np.append(new_str,atttime_start[i])
					new_stp=np.append(new_stp,sgti_stp[j])
					print 2,": gtistp in"
			else:
				if (atttime_stop[i]>=sgti_str[j])&(atttime_stop[i]<=sgti_stp[j]):
					new_str=np.append(new_str,sgti_str[j])
					new_stp=np.append(new_stp,atttime_stop[i])
					print 3,": gtistr in"
				else:
					new_str=np.append(new_str,sgti_str[j])
					new_stp=np.append(new_stp,sgti_stp[j])
					print 4,": gti in"
	
	print "gti sum: ",(new_stp-new_str).sum()
	if (new_stp-new_str).sum()>0:
		tbhdu=[[] for j in range(lenfile-2)]
		
		for i in range(1,lenfile-1):
			col1 = pf.Column(name = 'START', format='D',array = new_str)
			col2 = pf.Column(name = 'STOP', format='D',array = new_stp)
			cols = pf.ColDefs([col1,col2])
			tbhdu[i-1] = pf.BinTableHDU.from_columns(cols,name = 'GTI%s'%(i-1))
		
		prihdu = pf.PrimaryHDU()
		list1=[prihdu]
		for k in range(len(tbhdu)):
			list1=np.append(list1,tbhdu[k])
		
		list1=np.append(list1,softhdu)
		list1=list(list1)
		tbhdulist = pf.HDUList(list1)
		tbhdulist.writeto(Wbkg_gtifile, clobber=True)
		return new_str,new_stp

def hptime(Wpath, ObsID, dettb, scrtype, sel):
	OrgWpath = Wpath + '/Org'  ###
	Org_screenpath = "%s/%s" % (OrgWpath, 'Screen')  ###

	l=0
	fig = plt.figure(figsize = (16,9))
	gs = gridspec.GridSpec(1, 1,bottom=0.05,top=0.96,left=0.05,right=0.95,hspace=0.05)
	gs0 = gridspec.GridSpecFromSubplotSpec(3, 1,subplot_spec=gs[0], hspace=0.)
	
	ax=[plt.subplot(gs0[0]),plt.subplot(gs0[1],sharex=plt.subplot(gs0[0])),plt.subplot(gs0[2],sharex=plt.subplot(gs0[0]))]
	tball=[[]]*2
	for dtb in dettb:
		n, bins, patches=ax[l].hist(dtb,int(dtb.max()-dtb.min())/10)
		bint=(bins[:-1]+bins[1:])/2
		if sel==1:
			dmedian=np.median(n[n>0])
			#dmedian=620
			#dstd=np.std(n[n>0])
			dmean=np.mean(n[n>0])
			#plt.show()
			#plt.close()
			m=np.argsort(n[n>0])
			dstd=np.std(n[n>0][m][:-60])
			#dstd=110
			stl=dstd
		else:
			dmedian=np.median(n[n>0])
			#dmedian=620
			#dstd=np.std(n[n>0])
			dmean=np.mean(n[n>0])
			#plt.show()
			#plt.close()
			m=np.argsort(n[n>0])
			dstd=np.std(n[n>0][m][:-60])
			stl=5*dstd
		bint1=bint[n>(dmedian+stl)][1:][(bint[n>(dmedian+stl)][1:]-bint[n>(dmedian+stl)][:-1])>20]
		bint2=bint[n>(dmedian+stl)][:-1][(bint[n>(dmedian+stl)][1:]-bint[n>(dmedian+stl)][:-1])>20]
		if bint1.size>0:
			tb01=np.r_[bint[n>(dmedian+stl)][0],bint1]-10
			tb02=np.r_[bint2,bint[n>(dmedian+stl)][-1]]+10
			tb1=tb01[(tb02-tb01)>20]
			tb2=tb02[(tb02-tb01)>20]
			tball[0]=np.append(tball[0],tb1)
			tball[1]=np.append(tball[1],tb2)
		
		ax[l].scatter(bint,n,s=7,c='k',label='box %s'%l)
		ax[l].scatter(bint[n>(dmedian+stl)],n[n>(dmedian+stl)],c='r',s=10)
		ax[l].grid(True)
		ax[l].legend()
		collection = collections.BrokenBarHCollection.span_where(bint, ymin=-1000, ymax=10000, where=np.array(n>(dmedian+stl),dtype=int) > 0, alpha=0.3, color='orange')
		ax[l].add_collection(collection)
		l=l+1
	
	ax[0].set_title("%s" % (ObsID))
	ax[0].tick_params(labelbottom=False)
	ax[1].tick_params(labelbottom=False)
	ax[1].set_ylabel('counts/10s (1 box, good time)',fontsize=12)
	ax[2].set_xlabel('Time (s)',fontsize=12)
	plt.savefig('%s/%s_tct_%s.png' % (Org_screenpath, ObsID, scrtype), dpi = 300)
	if len(tball[0])>0:
		tbnew=np.array(tball).T[np.lexsort(tball[::-1])].T
		gtnew0=len(tbnew[0])
		jt1=tbnew[0][0]
		jt2=tbnew[1][0]
		gtnew=[[]]*2
		for i in range(1,gtnew0):
			if jt2<tbnew[0][i]:
				gtnew[0]=np.append(gtnew[0],jt1)
				gtnew[1]=np.append(gtnew[1],jt2)
				jt1=tbnew[0][i]
				jt2=tbnew[1][i]
			else:
				if jt2<tbnew[1][i]:
					jt2=tbnew[1][i]
		
		gtnew[0]=np.append(gtnew[0],jt1)
		gtnew[1]=np.append(gtnew[1],jt2)
		np.savetxt('%s/%s_hptime_%s.txt' % (Org_screenpath, ObsID, scrtype), gtnew, fmt="%s")
		return gtnew
	else:
		gtnew=0
		return gtnew

def read_scr(Wpath, ObsID, scrtype):
	OrgWpath = Wpath + '/Org'  ###
	Org_screenpath = "%s/%s" % (OrgWpath, 'Screen')  ###
	scdata=pf.open("%s/%s_me_screen_%s.fits" % (Org_screenpath, ObsID, scrtype))#me_screen_bkg.fits
	evtdata=scdata[1].data
	detidor=evtdata.field(1)
	#ch=evtdata.field(2)
	#etype=evtdata.field(6)
	#flag=evtdata.field(7)
	pior=evtdata.field(8)
	#grade=evtdata.field(9)
	gtidata=scdata[2].data
	startt=gtidata.field(0)
	stopt=gtidata.field(1)
	head_sc=scdata[1].header
	tstime=head_sc['TSTART']
	gtime=np.sum(stopt-startt)
	dett0=evtdata.field(0)
	dettb0=dett0[(pior<631)&(pior>68)&(detidor<576)]
	dettb1=dett0[(pior<631)&(pior>68)&(detidor>=576)&(detidor<1152)]
	dettb2=dett0[(pior<631)&(pior>68)&(detidor>=1152)]
	dettb=[dettb0,dettb1,dettb2]
	return detidor, pior, startt, stopt, tstime, gtime, dett0, dettb

def mescreen_newbd(Wpath, ObsID,Program_dir):
	OrgWpath = Wpath + '/Org'  ###
	Org_obspath = "%s/%s" % (OrgWpath, ObsID)  ###
	Org_gtipath = "%s/%s" % (OrgWpath, 'GTI')  ###
	Org_screenpath = "%s/%s" % (OrgWpath, 'Screen')  ###
	Wgradefile = "%s/%s_me_grade.fits" % (Org_obspath, ObsID)  ###
	Wbkg_gtifile = "%s/%s_me_gti_bkg.fits" % (Org_gtipath, ObsID)  ###
	Wbkg_screenfile = "%s/%s_me_screen_bkg.fits" % (Org_screenpath, ObsID)  ###
	Wbaddetfile = "%s/%s_me_bdet.fits" % (Org_obspath, ObsID)  ###bad det
	Wbaddeterrfile = "%s/%s_me_bderr.fits" % (Org_obspath, ObsID)  ###bad det
	Wrebin_screenfile = "%s/%s_me_screen_rebin.fits" % (Org_screenpath, ObsID)  ###% s / % s / me_screen_rebin.fits

	bdet0file = '/hxmt/home/saina/tar/bd0.fits'
	#Wpath='/hbkg/user/saina/data294/P0211007/'
	#ObsID='P0211007131'
	###Program_dir = '/sharefs/hbkg/user/saina/tar'###

	file_bdpath = '%s/%s'%(Program_dir,ObsID[:8])
	path= Wpath +'/' + ObsID###+"/"
	ft=25
###############################################pixel_class##############################################
	data = xlrd.open_workbook("%s/Blank_sky/pixel_class_mu.xlsx"%(Program_dir))
	classall=data.sheets()[3]
	class_p2=data.sheets()[0]
	class_st=data.sheets()[1]
	class_p3=data.sheets()[2]

	num_p2=class_p2.col_values(0)
	num_p2_e1=class_p2.col_values(1)
	num_p2_e2=class_p2.col_values(2)
	num_p2_e3=class_p2.col_values(3)
	num_p2_err=class_p2.col_values(4)
	num_p2_mu1=class_p2.col_values(6)
	num_p2_mu2=class_p2.col_values(7)

	num_st=class_st.col_values(0)
	num_st_cutp=class_st.col_values(1)
	num_st_cutt=class_st.col_values(2)
	num_st_e1=class_st.col_values(3)
	num_st_e2=class_st.col_values(4)
	num_st_e3=class_st.col_values(5)
	num_st_err=class_st.col_values(6)
	num_st_mu1=class_st.col_values(8)
	num_st_mu2=class_st.col_values(9)

	num_p3=class_p3.col_values(0)
	num_p3_e1=class_p3.col_values(1)
	num_p3_e2=class_p3.col_values(2)
	num_p3_e3=class_p3.col_values(3)
	num_p3_e4=class_p3.col_values(4)
	num_p3_err=class_p3.col_values(5)
	num_p3_mu1=class_p3.col_values(7)
	num_p3_mu2=class_p3.col_values(8)

	num=classall.col_values(0)
	classnum=classall.col_values(3)

	class1=np.array(num)[np.array(classnum)<=1]
	class2=np.array(num)[(np.array(classnum)<=2)&(np.array(classnum)>1)]
	class3=np.array(num)[(np.array(classnum)<=3)&(np.array(classnum)>2)]
	class4=np.array(num)[(np.array(classnum)<=4)&(np.array(classnum)>3)]
	class5=np.array(num)[(np.array(classnum)<=5)&(np.array(classnum)>4)]
	badclass=np.array(num)[(np.array(classnum)>5)]

###########################################high point/bkg gti####################################
	new_str,new_stp=bkg_gti(Wpath, ObsID, ft,Program_dir)

#####################################mescreen(all pixel)#########################################
	if (len(new_str)!=0)&((new_stp-new_str).sum()>300):
		#os.environ["PFILES"]="/sharefs/hbkg/user/saina/pfiles%s;/home/hxmt/hxmtsoft2/soft/install/x86_64-unknown-linux-gnu-libc2.12/syspfiles"%(ObsID)
		cmd = 'mescreen evtfile=%s gtifile=%s outfile=%s baddetfile=%s userdetid="0-53"' % (Wgradefile, Wbkg_gtifile, Wbkg_screenfile,bdet0file)
		print('mescreen(all pixel) start')
		print(cmd)
		os.system(cmd)
		print('mescreen(all pixel) end')
#####################################idcount & hr################################################
		scrtype='bkg'
		sel=1
		detidor, pior, startt, stopt, tstime, gtime, dett0, dettb = read_scr(Wpath, ObsID, scrtype)
		gtnew=hptime(Wpath, ObsID, dettb, scrtype, sel)
		print("gtnew:", gtnew)
		if gtnew==0:
			detid=detidor
			pi=pior
		else:
			condi1 =  ((new_stp-new_str).sum()-(gtnew[1]-gtnew[0]).sum())
			print("(new_stp-new_str).sum()-(gtnew[1]-gtnew[0]).sum():",condi1)
			if condi1<=250:
				file_bd = '%s/run_merun.txt'%(file_bdpath)###
				with open(file_bd,'a') as f:
					f.write('%s -- Warning: too short gti -- %s' % (ObsID, ((new_stp - new_str).sum() - (gtnew[1] - gtnew[0]).sum())) + '\n')
				for i in range(int(ObsID[7:])):
					###obid1='P021100'+str(int(ObsID[7:]) - i)
					obid1 = 'P030124062901'  ###change by LQ
					Org_obid1path = "%s/%s" % (OrgWpath, obid1)
					if os.path.exists("%s/%s_me_bdet.fits" % (Org_obid1path,obid1)):###careful
						os.system('cp %s/%s_me_bdet.fits %s' % (Org_obid1path,obid1, Wbaddetfile))###
						with open(file_bd,'a') as f:
							f.write('%s -- cp bdfile from %s' % (ObsID, obid1) + '\n')
						print "cp bdfile from ",obid1
						break

				return 0

			mm1=(dett0<gtnew[0][0])|(dett0>gtnew[1][0])
			if len(gtnew[0])>2:
				for i in range(1,len(gtnew[0])):
					mm2=mm1&(dett0<gtnew[0][i])|(dett0>gtnew[1][i])
					mm1=mm2
			else:
				mm2=mm1

			detid=detidor[mm2]
			pi=pior[mm2]
		print("idcount & hr")
		print("detid.size:", detid.size)
###########################################################################################
		idcount_t=[]
		idcount=[]
		hd1=[]
		hd11=[]
		hd12=[]
		for i in range(1728):
			detidx=detid[(detid==i)&(pi<631)&(pi>68)]
			idcount_t=np.append(idcount_t,len(detidx)/gtime)
			idcount=np.append(idcount,len(detidx))
			detid1x1=detid[(detid==i)&(pi<102)&(pi>68)]
			detid1x2=detid[(detid==i)&(pi<137)&(pi>102)]
			hd11=np.append(hd11,len(detid1x1))
			hd12=np.append(hd12,len(detid1x2))
			if (len(detid1x1)*len(detid1x2))!=0:
				hd1=np.append(hd1,float(len(detid1x1))/float(len(detid1x2)))
			else:
				hd1=np.append(hd1,0)
		print("idcount_t,hd11,hd12:", idcount_t,hd11,hd12)
#########################################ideal ct##############################################

		sbd=pf.open("/home/hxmt/hxmtsoft2/soft/install/x86_64-unknown-linux-gnu-libc2.12/refdata/medetectorstatus.fits")
		nump=range(1728)
		#####################################1: p2 ict(only 124)########
		#class12=np.r_[class1,class2]
		class124=np.r_[class1,class2,class4]

		maskp2=np.in1d(nump,class124)
		mask=np.in1d(num_p2,class124)

		num1=np.array(num_p2)[mask]
		el1=np.array(num_p2_e1)[mask]
		el2=np.array(num_p2_e2)[mask]
		el3=np.array(num_p2_e3)[mask]
		er=np.array(num_p2_err)[mask]
		mu1=np.array(num_p2_mu1)[mask]
		mu2=np.array(num_p2_mu2)[mask]

		txtime=[tstime]*len(num1)
		tx=(txtime-mu1)/mu2

		ictup_p2=el1*(tx**2)+el2*tx+el3+np.abs(er)*3
		ict_p2=el1*(tx**2)+el2*tx+el3

		##################################2: st ict#######
		maskst=np.in1d(nump,class3)
		part1=[[]]*len(class3)
		allcut=[]*len(class3)
		for i in range(len(class3)):
			part1[i]=class_st.row_values(i*2)
			allcut=np.append(allcut,class_st.row_values(i*2)[2])

		part2=[[]]*len(class3)
		for i in range(len(class3)):
			part2[i]=class_st.row_values(i*2+1)

		txs=el31=el32=el33=err3=[]
		for i in range(len(class3)):
			if tstime<allcut[i]:
				txs=np.append(txs,(tstime-part1[i][8])/part1[i][9])
				el31=np.append(el31,part1[i][3])
				el32=np.append(el32,part1[i][4])
				el33=np.append(el33,part1[i][5])
				err3=np.append(err3,part1[i][6])
			else:
				txs=np.append(txs,(tstime-part2[i][8])/part2[i][9])
				el31=np.append(el31,part2[i][3])
				el32=np.append(el32,part2[i][4])
				el33=np.append(el33,part2[i][5])
				err3=np.append(err3,part2[i][6])

		ictup_st=el31*(txs**2)+el32*txs+el33+np.abs(err3)*3
		ict_st=el31*(txs**2)+el32*txs+el33

		##################################3: p3 ict(if)###########
		maskp3=np.in1d(nump,class5)
		mask3=np.in1d(num_p3,class5)

		p3time=[tstime]*len(class5)
		tx3=(p3time-np.array(num_p3_mu1)[mask3])/np.array(num_p3_mu2)[mask3]

		ict_p3=np.array(num_p3_e1)[mask3]*(tx3**3)+np.array(num_p3_e2)[mask3]*(tx3**2)+np.array(num_p3_e4)[mask3]+np.array(num_p3_e3)[mask3]*tx3
		ictup_p3=ict_p3+np.abs(np.array(num_p3_err)[mask3])*3

		##################################comb######################
		maskp=maskp2|maskst
		num_all=np.r_[num1,class3]
		ictup_all=np.r_[ictup_p2,ictup_st]
		ict_all=np.r_[ict_p2,ict_st]
		newl=np.array([list(num_all),list(ictup_all),list(ict_all)])
		newl1=newl.T[np.lexsort(newl[::-1,:])].T
		#####################################pic####################
		#maskp124=np.in1d(nump,class124)

		'''plt.plot(nump,idcount_t)
		plt.scatter(newl1[0],newl1[1],s=1,c='r')
		plt.scatter(newl1[0][(newl1[1]<idcount_t[maskp])|(idcount_t[maskp]>0.1)],newl1[1][(newl1[1]<idcount_t[maskp])|(idcount_t[maskp]>0.1)],s=5,c='k')
		plt.scatter(newl1[0][hd1[maskp]>3],newl1[1][hd1[maskp]>3],s=5,c='g')
		plt.ylim(0,0.2)
		plt.show()'''
#####################################compare ict&idcount###############################

		selp1=((newl1[1]>=idcount_t[maskp])&(idcount_t[maskp]<=0.1))
		#selp2=(newl1[1]<idcount_t[maskp])&(idcount_t[maskp]<newl1[2]+0.02)&(np.in1d(newl1[0],class1))
		selp3=(hd1[maskp]<=3)
		#selp=(selp1|selp2)&selp3
		selp=selp1&selp3
		bdm=np.in1d(np.array(nump),newl1[0][selp],invert=True)
		newbd=np.array(nump)[bdm]
		numbd=len(newbd)
		print("newbd,numbd:", newbd,numbd)
		#plt.scatter(np.array(nump),idcount_t)
		#plt.scatter(np.array(nump)[maskp],ictup)
################################bad detector list(ObsID)###################################

		head_sbd=sbd[1].header

		col1 = pf.Column(name = 'DetID', format='1I',unit='Det ID',array = newbd)
		col2 = pf.Column(name = 'TIMERANGE', format='20A',unit='time',array = np.array([0]*numbd))
		col3 = pf.Column(name = 'TIMERANGE2', format='20A',unit='time',array = np.array(['INDEF']*numbd))
		col4 = pf.Column(name = 'TYPE', format='20A',unit='bad/hot/flick/warm/cold',array = np.array(['Bad']*numbd))
		col5 = pf.Column(name = 'STATUS', format='1B',unit='bad',array = np.array([0]*numbd))
		cols1 = pf.ColDefs([col1,col2,col3,col4,col5])
		tbhdu0 = pf.BinTableHDU.from_columns(cols1,name = 'detectorStatus')
		tbhdu0.header=head_sbd
		prihdu = pf.PrimaryHDU()
		tbhdulist = pf.HDUList([prihdu,tbhdu0])
		tbhdulist.writeto(Wbaddetfile, clobber=True)###
		file_bd = '%s/run_merun.txt' % file_bdpath
		plt.close('all')
		if numbd>480:
			with open(file_bd,'a') as f:
				f.write('%s -- Warning: too much bad detector -- %s' % (ObsID, numbd) + '\n')
			print 'Warning: too much bad detector -- ',numbd
			if numbd>580:
				os.system("mv %s %s" % (Wbaddetfile, Wbaddeterrfile))###
				print('i range',ObsID[7:8])
				obid1 = 'P030124062901'###change by LQ
				Org_obid1path = "%s/%s" % (OrgWpath, obid1)
				os.system('cp %s/%s_me_bdet.fits %s' % (Org_obid1path, obid1, Wbaddetfile))
				with open(file_bd, 'a') as f:
					f.write('%s -- cp bdfile from %s' % (ObsID, obid1) + '\n')
				print "cp bdfile from ", obid1
				'''###I don't understand the logic of follow
				for i in range(int(ObsID[7:8])):
					obid1= ObsID[:7] + str(int(ObsID[7:11]) - i) + '01'
					if os.path.exists('%s/%s_me_bdet.fits'%(Org_obspath, obid1)):
						os.system('cp %s/%s_me_bdet.fits %s' % (Org_obspath, obid1, Wbaddetfile))
						with open(file_bd,'a') as f:
							f.write('%s -- cp bdfile from %s' % (ObsID, obid1) + '\n')
						print "cp bdfile from ",obid1
						break
				'''
	else:
		file_bd = '%s/run_merun.txt'%(file_bdpath)
		with open(file_bd,'a') as f:
				f.write('%s -- Warning: too short gti -- %s' % (ObsID, (new_stp - new_str).sum()) + '\n')
		for i in range(int(ObsID[7:])):
			obid1= ObsID[:7] + str(int(ObsID[7:11]) - i) + '01'
			if os.path.exists('%s/%s_me_bdet.fits'%(Org_obspath, obid1)):
				os.system('cp %s/%s_me_bdet.fits %s' % (Org_obspath, obid1, Wbaddetfile))
				with open(file_bd,'a') as f:
					f.write('%s -- cp bdfile from %s' % (ObsID, obid1) + '\n')
				print "cp bdfile from ",obid1
				break

###########################################hp2(in gen lc)##################################################
	scrtype2='rebin'
	sel2=2
	os.system('mescreen evtfile=%s \
		gtifile=%s outfile=%s \
		baddetfile=%s userdetid="0-53"' % (Wgradefile, Wbkg_gtifile, Wrebin_screenfile, Wbaddetfile))
	detidor, pior, startt, stopt, tstime, gtime, dett0, dettbr = read_scr(Wpath, ObsID, scrtype2)
	gtnew_re=hptime(Wpath, ObsID, dettbr, scrtype2, sel2)

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
	###path = '/hxmt/work/HXMT-DATA/1L/A%s/%s/' % (ObsID[1:3], ObsID[:-5])

	###Wpath='/sharefs/hbkg/user/saina/data294/P0201013/'
	#lclist= os.listdir(Wpath)
	#lclist.sort()
	###lclist='P020101359901'
	#Wpath=Wpath+lclist
	#ehkfile=Wpath+lclist+"/ehk.fits"
	#gtifile=Wpath+'me_gti.fits'
	mescreen_newbd(Wpath,ObsID,program_tree)
