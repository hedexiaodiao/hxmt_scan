#!/hxmt/soft/Develop/anaconda2/bin/python
import numpy as np
import astropy.io.fits as pf
import math
import os
import time
from collections import Counter
from astropy import units as u
from astropy.coordinates import SkyCoord
import glob

def mkdir_try(dirname):
    if os.path.exists(dirname) == 0:
        try:
            os.makedirs(dirname)
        except OSError:
            print('Wrong with make dir:\n'+dirname)

def dataout(Wpath, ObsID):
	OrgWpath = Wpath + '/Org'  ###
	NetWpath = Wpath + '/Net'  ###
	Net_obspath = "%s/%s" % (NetWpath, ObsID)  ###
	Org_attpath = "%s/%s" % (OrgWpath, 'Att')  ###
	Wattfile = "%s/%s_Att.fits" % (Org_attpath, ObsID)  ###
	Midd_obspath = Wpath + '/Midd/' + ObsID  ###
	fProd_obspath = Wpath + '/Prod/' + ObsID  ###
	path = Wpath + ObsID + "/ME/"
	if not os.path.exists("%s/%s_ub_chi2.dat"%(Midd_obspath, ObsID)):
		print "%s: fit not end"%Midd_obspath
		#continue
		return 0
	else:
		print "%s: process start"%(Midd_obspath)
		################################source1#######################################
		resdata = np.loadtxt("%s/%s_data_3box.dat" % (Midd_obspath, ObsID))
		bkgfile = np.loadtxt("%s/%s_chi2.dat" % (Midd_obspath, ObsID))
		res0 = resdata[:,3]
		res1 = resdata[:,7]
		res2 = resdata[:,11]
		model0 = resdata[:,2]
		model1 = resdata[:,6]
		model2 = resdata[:,10]
		nolen = len(res0)
		print "res len: ",len(res0)
		st0=bkgfile[3]
		#count0=Counter(model0).most_common(1)
		#st0=count0[0][0]
		#count1=Counter(model1).most_common(1)
		#st1=count1[0][0]
		#count2=Counter(model2).most_common(1)
		#st2=count2[0][0]
		#bkg_st = [st0,st1,st2]
		bkg_st=[bkgfile[3],bkgfile[6],bkgfile[9]]
		bkg_name = ['BACKGRND_box1','BACKGRND_box2','BACKGRND_box3']
		f_bkg = [None,None,None]
		bkg_num = [0.,0.,0.]
		bkg_num1 = [bkgfile[4],bkgfile[7],bkgfile[10]]
		bkg_num2 = [bkgfile[5],bkgfile[8],bkgfile[11]]
		if os.path.exists("%s/%s_source_moreinfo.txt"%(Midd_obspath, ObsID)):
			src_list_all = np.loadtxt("%s/%s_source_moreinfo.txt" % (Midd_obspath, ObsID), delimiter='\t', dtype=bytes).astype(str)
			if len(src_list_all.shape)>1:
				src_err1_all = src_list_all[:,6].astype(float)
				src_list = src_list_all[src_list_all[:,6].astype(float)>0]
			else:
				src_err1_all = [src_list_all[6].astype(float)]
				src_list = src_list_all
			#print src_list
			if len(src_list.shape)>1:
				#src_num1 = src_list[:,0].astype(int)
				src_num1 = range(1,len(src_list)+1)
				src_num2 = src_list[:,1].astype(int)
				src_name = src_list[:,2]
				src_ra = src_list[:,3].astype(float)
				src_dec = src_list[:,4].astype(float)
				src_nor = src_list[:,5].astype(float)
				#src_err = src_list[:,6].astype(float)
				#src_snr = src_list[:,7].astype(float)
				src_err1 = src_list[:,6].astype(float)
				src_err2 = src_list[:,7].astype(float)
				src_snr = src_list[:,9].astype(float)
			else:
				src_num1_all=1
				if src_err1_all[0]>0:
					src_num1 = [src_list[0].astype(int)]
					src_num2 = [src_list[1].astype(int)]
					src_name = [src_list[2]]
					src_ra = [src_list[3].astype(float)]
					src_dec = [src_list[4].astype(float)]
					src_nor = [src_list[5].astype(float)]
					#src_err = [src_list[6].astype(float)]
					#src_snr = [src_list[7].astype(float)]
					src_err1 = [src_list[6].astype(float)]
					src_err2 = [src_list[7].astype(float)]
					src_snr = [src_list[9].astype(float)]
				else:
					src_num1=src_num2=src_name=src_ra=src_dec=src_nor=src_err1=src_err2=src_snr=[]
		else:
				src_num1=src_num2=src_name=src_ra=src_dec=src_nor=src_err1=src_err2=src_snr=[]
				src_err1_all=[]
		################################source2#######################################
		if os.path.exists("%s/%s_new.txt"%(Midd_obspath, ObsID)):
			usrc_list = np.loadtxt("%s/%s_new.txt" % (Midd_obspath, ObsID))
			if len(usrc_list.shape)>1:
				usrc_num1 = [usrc_list[:,0][-1]]
				usrc_num2 = [([0]*len(usrc_list.shape))[-1]]
				usrc_ra = [usrc_list[:,1][-1]]
				usrc_dec = [usrc_list[:,3][-1]]
				usrc_nor = [usrc_list[:,5][-1]]
				#src_err = src_list[:,6]
				#src_snr = src_list[:,7]
				usrc_err1 = [usrc_list[:,6][-1]]
				usrc_err2 = [usrc_list[:,7][-1]]
				usrc_snr = [usrc_list[:,9][-1]]
				usrc_name =[]
				for i in range(len(usrc_num1)):
					us_icrs = SkyCoord(ra=usrc_ra[i]*u.degree, dec=usrc_dec[i]*u.degree, frame='icrs')
					nr1=abs(int(us_icrs.ra.hms.h))
					nr2=abs(int(us_icrs.ra.hms.m))
					if nr2<10:
						nr2='0'+str(nr2)
					if nr1<10:
						nr1='0'+str(nr1)
					np1='ISTC J'+str(nr1)+str(nr2)
					nd1=abs(int(us_icrs.dec.hms.h))
					nd2=abs(int(us_icrs.dec.hms.m))
					if nd2<10:
						nd2='0'+str(nd2)
					if nd1<10:
						nd1='0'+str(nd1)
					if int(us_icrs.dec.hms.h)<0:
						np2='-'+str(nd1)+str(nd2)
					else:
						np2='+'+str(nd1)+str(nd2)
					usrc_name = np.append(usrc_name,np1+np2)
			else:
				usrc_num1 = [usrc_list[0]]
				usrc_num2 = [0]
				#usrc_name = [None]
				usrc_ra = [usrc_list[1]]
				usrc_dec = [usrc_list[3]]
				usrc_nor = [usrc_list[5]]
				#src_err = [src_list[6]]
				#src_snr = [src_list[7]]
				usrc_err1 = [usrc_list[6]]
				usrc_err2 = [usrc_list[7]]
				usrc_snr = [usrc_list[9]]
				us_icrs = SkyCoord(ra=usrc_ra*u.degree, dec=usrc_dec*u.degree, frame='icrs')
				nr1=abs(int(us_icrs.ra.hms.h))
				nr2=abs(int(us_icrs.ra.hms.m))
				if nr2<10:
					nr2='0'+str(nr2)
				if nr1<10:
					nr1='0'+str(nr1)
				np1='ISTC J'+str(nr1)+str(nr2)
				nd1=abs(int(us_icrs.dec.hms.h))
				nd2=abs(int(us_icrs.dec.hms.m))
				if nd2<10:
					nd2='0'+str(nd2)
				if nd1<10:
					nd1='0'+str(nd1)
				if int(us_icrs.dec.hms.h)<0:
					np2='-'+str(nd1)+str(nd2)
				else:
					np2='+'+str(nd1)+str(nd2)
				usrc_name = [np1+np2]
		else:
				usrc_num1=usrc_num2=usrc_name=usrc_ra=usrc_dec=usrc_nor=usrc_err1=usrc_err2=usrc_snr=[]
		################################source3#######################################
		if os.path.exists("%s/%s_ub_moreinfo.txt"%(Midd_obspath, ObsID)):
			wsrc_list = np.loadtxt("%s/%s_ub_moreinfo.txt" % (Midd_obspath, ObsID), delimiter='\t', dtype=bytes).astype(str)
			i2=0
			for i1 in range(len(wsrc_list)):
				wsrc_list[i2][0]=i2+1
				if (wsrc_list[i2][9].astype(float)==99999999)|(wsrc_list[i2][9]==None)|\
				((wsrc_list[i2][9].astype(float)==99999999)&(wsrc_list[i2][5].astype(float)<=0))|\
				((wsrc_list[i2][9].astype(float)<=0)&(wsrc_list[i2][5].astype(float)<=0))|\
				((wsrc_list[i2][5].astype(float)<=-22)):
					wsrc_list=np.delete(wsrc_list,i2,axis=0)
					i2=i2-1
				i2=i2+1
			if len(wsrc_list.shape)>1:
				model_num=wsrc_list[:,0].astype(int)
				wsrc_num1 = range(1,len(wsrc_list)+1)
				wsrc_num2 = wsrc_list[:,1].astype(int)
				wsrc_name = wsrc_list[:,2]
				wsrc_ra = wsrc_list[:,3].astype(float)
				wsrc_dec = wsrc_list[:,4].astype(float)
				wsrc_nor = wsrc_list[:,5].astype(float)
				#src_err = src_list[:,6].astype(float)
				#src_snr = src_list[:,7].astype(float)
				wsrc_err1 = wsrc_list[:,6].astype(float)
				wsrc_err2 = wsrc_list[:,7].astype(float)
				wsrc_snr = wsrc_list[:,9].astype(float)
			else:
				model_num= [wsrc_list[0].astype(int)]
				wsrc_num1 = [int(1)]
				wsrc_num2 = [wsrc_list[1].astype(int)]
				wsrc_name = wsrc_list[2]
				wsrc_ra = [wsrc_list[3].astype(float)]
				wsrc_dec = [wsrc_list[4].astype(float)]
				wsrc_nor = [wsrc_list[5].astype(float)]
				#src_err = [src_list[6].astype(float)]
				#src_snr = [src_list[7].astype(float)]
				wsrc_err1 = [wsrc_list[6].astype(float)]
				wsrc_err2 = [wsrc_list[7].astype(float)]
				wsrc_snr = [wsrc_list[9].astype(float)]
		else:
				wsrc_num1=wsrc_num2=wsrc_name=wsrc_ra=wsrc_dec=wsrc_nor=wsrc_err1=wsrc_err2=wsrc_snr=[]
		b1_name = bkg_name+list(src_name)+list(usrc_name)+list(wsrc_name)
		b1_num1 = bkg_num+list(src_num1)+list(usrc_num1)+list(wsrc_num1)
		b1_num2 = bkg_num+list(src_num2)+usrc_num2+list(wsrc_num2)
		b1_ra = f_bkg+list(src_ra)+list(usrc_ra)+list(wsrc_ra)
		b1_dec = f_bkg+list(src_dec)+list(usrc_dec)+list(wsrc_dec)
		b1_nor = bkg_st+list(src_nor)+list(usrc_nor)+list(wsrc_nor)
		b1_noru = bkg_num1+list(src_err1)+list(usrc_err1)+list(wsrc_err1)
		b1_norl = bkg_num2+list(src_err2)+list(usrc_err2)+list(wsrc_err2)
		#b1_noru = f_bkg+list(src_err1)
		#b1_norl = f_bkg+list(src_err2)
		b1_snr = f_bkg+list(src_snr)+list(usrc_snr)+list(wsrc_snr)
		col1 = pf.Column(name = 'NAME', format='20A',array = b1_name)
		col2 = pf.Column(name = 'FITORDER', format='I',array = b1_num1)
		col3 = pf.Column(name = 'SEQUENSE', format='I',array = b1_num2)
		col4 = pf.Column(name = 'RA', format='E',array = b1_ra)
		col5 = pf.Column(name = 'DEC', format='E',array = b1_dec)
		col6 = pf.Column(name = 'NORM', format='E',array = b1_nor)
		col7 = pf.Column(name = 'NORM_UB', format='E',array = b1_noru)
		col8 = pf.Column(name = 'NORM_LB', format='E',array = b1_norl)
		col9 = pf.Column(name = 'SNR', format='E',array = b1_snr)
		cols1 = pf.ColDefs([col1,col2,col3,col4,col5,col6,col7,col8,col9])
		tbhdu1 = pf.BinTableHDU.from_columns(cols1,name='SOURCES')
		######################################lc*3###########################################
		slcfile_list=glob.glob(r'%s/*small_g0*.lc'%(Net_obspath))###Warning, maybe wrong
		#blcfile_list=glob.glob(r'%s/*blind_g0*.lc'%(Wpath))
		slcfile_list.sort()
		#blcfile_list.sort()
		lc_o0 = pf.open(slcfile_list[0])
		lc_o1 = pf.open(slcfile_list[1])
		lc_o2 = pf.open(slcfile_list[2])
		o0_data = lc_o0[1].data
		o1_data = lc_o1[1].data
		o2_data = lc_o2[1].data
		o0t=o0_data.field(0)
		o0lc=o0_data.field(1)
		o1t=o1_data.field(0)
		o1lc=o1_data.field(1)
		o2t=o2_data.field(0)
		o2lc=o2_data.field(1)
		lc_all0 = pf.open("%s/me_lc_box0_small_cut.fits"%(Net_obspath))
		lc_all1 = pf.open("%s/me_lc_box1_small_cut.fits"%(Net_obspath))
		lc_all2 = pf.open("%s/me_lc_box2_small_cut.fits"%(Net_obspath))
		lc_data0 = lc_all0[1].data
		lc_gti = lc_all0[2]
		lc_time0 = lc_data0['Time']
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
		mask0=np.in1d(o0t,lc_time0)
		mask1=np.in1d(o1t,lc_time1)
		mask2=np.in1d(o2t,lc_time2)
		o0lcn=o0lc[mask0]
		o1lcn=o1lc[mask1]
		o2lcn=o2lc[mask2]
		bkg0=o0lcn-lc_ct0
		bkg1=o1lcn-lc_ct1
		bkg2=o2lcn-lc_ct2
		stb20=[bkgfile[3]]*nolen
		stb21=[bkgfile[6]]*nolen
		stb22=[bkgfile[9]]*nolen
		col01 = pf.Column(name = 'TIME', format='1D',array = lc_time0)
		col02 = pf.Column(name = 'ORG_CTS', format='E',array = o0lcn)
		col03 = pf.Column(name = 'SKY_BKG', format='E',array = bkg0)
		col04 = pf.Column(name = 'NET_CTS', format='E',array = lc_ct0)
		col05 = pf.Column(name = 'ORG_CTS_ERR', format='E',array = lc_ct_e0)
		col06 = pf.Column(name = 'PSF_SUM', format='E',array = model0)
		col07 = pf.Column(name = 'RESIDUAL', format='E',array = res0)
		col08 = pf.Column(name = 'BKG_CONSTANT', format='E',array = stb20)
		
		col11 = pf.Column(name = 'TIME', format='1D',array = lc_time1)
		col12 = pf.Column(name = 'ORG_CTS', format='E',array = o1lcn)
		col13 = pf.Column(name = 'SKY_BKG', format='E',array = bkg1)
		col14 = pf.Column(name = 'NET_CTS', format='E',array = lc_ct1)
		col15 = pf.Column(name = 'ORG_CTS_ERR', format='E',array = lc_ct_e1)
		col16 = pf.Column(name = 'PSF_SUM', format='E',array = model1)
		col17 = pf.Column(name = 'RESIDUAL', format='E',array = res1)
		col18 = pf.Column(name = 'BKG_CONSTANT', format='E',array = stb21)
		
		col21 = pf.Column(name = 'TIME', format='1D',array = lc_time2)
		col22 = pf.Column(name = 'ORG_CTS', format='E',array = o2lcn)
		col23 = pf.Column(name = 'SKY_BKG', format='E',array = bkg2)
		col24 = pf.Column(name = 'NET_CTS', format='E',array = lc_ct2)
		col25 = pf.Column(name = 'ORG_CTS_ERR', format='E',array = lc_ct_e2)
		col26 = pf.Column(name = 'PSF_SUM', format='E',array = model2)
		col27 = pf.Column(name = 'RESIDUAL', format='E',array = res2)
		col28 = pf.Column(name = 'BKG_CONSTANT', format='E',array = stb22)
		ppsf0=[[]]*(len(src_num1)+len(usrc_num1)+len(wsrc_num1))
		ppsf1=[[]]*(len(src_num1)+len(usrc_num1)+len(wsrc_num1))
		ppsf2=[[]]*(len(src_num1)+len(usrc_num1)+len(wsrc_num1))
		colsrc0=[[]]*(len(src_num1)+len(usrc_num1)+len(wsrc_num1))
		colsrc1=[[]]*(len(src_num1)+len(usrc_num1)+len(wsrc_num1))
		colsrc2=[[]]*(len(src_num1)+len(usrc_num1)+len(wsrc_num1))
		#print "res:",len(res0),"lc:",len(lc_time0)
		if len(res0)!=len(lc_time0):
			print "Wrong files, please check the time of lc files and model files."
			return 0
			#continue
		else:
			flag=0
			modnum=0
			modnum_w=0
			k1=0
			k2=0
			i1=0
			print "all mod: ",len(src_num1)+len(usrc_num1)+len(wsrc_num1)
			for j in range(len(src_err1_all)):
				#time.sleep(5)
				modnum=modnum+1
				if src_err1_all[j]<=0:
					continue
				#else:
				#allm0 = resdata[:,2]
				#allm1 = resdata[:,6]
				#allm2 = resdata[:,10]
				flag=1
				mod=np.loadtxt("%s/model_%s.dat"%(Midd_obspath,modnum))
				mod0=mod[:,2]
				mod1=mod[:,5]
				mod2=mod[:,8]
				ppsf0[i1]=mod0
				ppsf1[i1]=mod1
				ppsf2[i1]=mod2
				colsrc0[i1]=pf.Column(name = '%s'%src_name[i1], format='E',array = ppsf0[i1])
				colsrc1[i1]=pf.Column(name = '%s'%src_name[i1], format='E',array = ppsf1[i1])
				colsrc2[i1]=pf.Column(name = '%s'%src_name[i1], format='E',array = ppsf2[i1])
				i1=i1+1
			print "src ObsID: ",modnum,"    src>0 ObsID: ",i1
			for i in range(i1,i1+len(usrc_num1)):
				flag=1
				modnum=modnum+1
				mod=np.loadtxt("%s/ns_model_%s.dat"%(Midd_obspath,modnum))
				mod0=mod[:,2]
				mod1=mod[:,5]
				mod2=mod[:,8]
				ppsf0[i]=mod0
				ppsf1[i]=mod1
				ppsf2[i]=mod2
				colsrc0[i]=pf.Column(name = '%s'%usrc_name[k1], format='E',array = ppsf0[i])
				colsrc1[i]=pf.Column(name = '%s'%usrc_name[k1], format='E',array = ppsf1[i])
				colsrc2[i]=pf.Column(name = '%s'%usrc_name[k1], format='E',array = ppsf2[i])
				k1=k1+1
			print "mod ObsID: ",modnum,"    us ObsID: ",k1
			for i in range(i1+len(usrc_num1),i1+len(usrc_num1)+len(wsrc_num1)):
				flag=1
				mod=np.loadtxt("%s/ub_model_%s.dat"%(Midd_obspath,model_num[modnum_w]))
				mod0=mod[:,2]
				mod1=mod[:,6]
				mod2=mod[:,10]
				ppsf0[i]=mod0
				ppsf1[i]=mod1
				ppsf2[i]=mod2
				colsrc0[i]=pf.Column(name = '%s'%wsrc_name[k2], format='E',array = ppsf0[i])
				colsrc1[i]=pf.Column(name = '%s'%wsrc_name[k2], format='E',array = ppsf1[i])
				colsrc2[i]=pf.Column(name = '%s'%wsrc_name[k2], format='E',array = ppsf2[i])
				k2=k2+1
				modnum_w=modnum_w+1
			print "ubs ObsID: ",modnum_w
				
			if flag==1:
				cols2 = pf.ColDefs([col01,col02,col03,col04,col05,col06,col07,col08]+list(colsrc0))
				tbhdu2 = pf.BinTableHDU.from_columns(cols2,name="DATA_BOX1")
				cols3 = pf.ColDefs([col11,col12,col13,col14,col15,col16,col17,col18]+list(colsrc1))
				tbhdu3 = pf.BinTableHDU.from_columns(cols3,name="DATA_BOX2")
				cols4 = pf.ColDefs([col21,col22,col23,col24,col25,col26,col27,col28]+list(colsrc2))
				tbhdu4 = pf.BinTableHDU.from_columns(cols4,name="DATA_BOX3")
			else:
				cols2 = pf.ColDefs([col01,col02,col03,col04,col05,col06,col07,col08])
				tbhdu2 = pf.BinTableHDU.from_columns(cols2,name="DATA_BOX1")
				cols3 = pf.ColDefs([col11,col12,col13,col14,col15,col16,col17,col18])
				tbhdu3 = pf.BinTableHDU.from_columns(cols3,name="DATA_BOX2")
				cols4 = pf.ColDefs([col21,col22,col23,col24,col25,col26,col27,col28])
				tbhdu4 = pf.BinTableHDU.from_columns(cols4,name="DATA_BOX3")
			
		##########################################att#########################################
		att_file=pf.open("%s"%Wattfile)
		attr=att_file[1].data
		attq=att_file[3].data
		att_time=attr.field(0)
		att_ra=attr.field(1)
		att_dec=attr.field(2)
		att_q1=attq.field(1)
		att_q2=attq.field(2)
		att_q3=attq.field(3)
		col51 = pf.Column(name = 'TIME', format='1D',array = att_time)
		col52 = pf.Column(name = 'RA', format='E',array = att_ra)
		col53 = pf.Column(name = 'DEC', format='E',array = att_dec)
		col54 = pf.Column(name = 'Q1', format='E',array = att_q1)
		col55 = pf.Column(name = 'Q2', format='E',array = att_q2)
		col56 = pf.Column(name = 'Q3', format='E',array = att_q3)
		cols5 = pf.ColDefs([col51,col52,col53,col54,col55,col56])
		tbhdu5 = pf.BinTableHDU.from_columns(cols5,name="ATT")
		
		#####################################blind################################################
		#blind_file0=pf.open("%s/me_0_blind_g0_10.lc"%Wpath)
		#blind_file1=pf.open("%s/me_1_blind_g0_28.lc"%Wpath)
		#blind_file2=pf.open("%s/me_2_blind_g0_46.lc"%Wpath)
		blcfile_list=glob.glob(r'%s/*blind_g0*.lc'%(Net_obspath))
		blcfile_list.sort()
		blind_file0=pf.open(blcfile_list[0])
		blind_file1=pf.open(blcfile_list[1])
		blind_file2=pf.open(blcfile_list[2])
		blind0=blind_file0[1].data
		blind1=blind_file1[1].data
		blind2=blind_file2[1].data
		mask0=np.in1d(blind0.field(0),lc_time0)
		mask1=np.in1d(blind1.field(0),lc_time1)
		mask2=np.in1d(blind2.field(0),lc_time2)
		blind0lc=blind0.field(1)[mask0]
		blind1lc=blind1.field(1)[mask1]
		blind2lc=blind2.field(1)[mask2]
		blind0exp=blind0.field(2)[mask0]
		blind1exp=blind1.field(2)[mask1]
		blind2exp=blind2.field(2)[mask2]
		blind0err=blind0.field(3)[mask0]
		blind1err=blind1.field(3)[mask1]
		blind2err=blind2.field(3)[mask2]
		blind0de=blind0.field(4)[mask0]
		blind1de=blind1.field(4)[mask1]
		blind2de=blind2.field(4)[mask2]
		col61 = pf.Column(name = 'TIME', format='1D',unit='s',array = lc_time0)
		col62 = pf.Column(name = 'COUNTS', format='E',unit='counts',array = blind0lc)
		col63 = pf.Column(name = 'FRACEXP', format='E',array = blind0exp)
		col64 = pf.Column(name = 'ERROR', format='E',unit='counts',array = blind0err)
		col65 = pf.Column(name = 'DEADC', format='E',array = blind0de)
		cols6 = pf.ColDefs([col61,col62,col63,col64,col65])
		tbhdu6 = pf.BinTableHDU.from_columns(cols6,name="BLINDDETECTOR1")
		col71 = pf.Column(name = 'TIME', format='1D',unit='s',array = lc_time1)
		col72 = pf.Column(name = 'COUNTS', format='E',unit='counts',array = blind1lc)
		col73 = pf.Column(name = 'FRACEXP', format='E',array = blind1exp)
		col74 = pf.Column(name = 'ERROR', format='E',unit='counts',array = blind1err)
		col75 = pf.Column(name = 'DEADC', format='E',array = blind1de)
		cols7 = pf.ColDefs([col71,col72,col73,col74,col75])
		tbhdu7 = pf.BinTableHDU.from_columns(cols7,name="BLINDDETECTOR2")
		col81 = pf.Column(name = 'TIME', format='1D',unit='s',array = lc_time2)
		col82 = pf.Column(name = 'COUNTS', format='E',unit='counts',array = blind2lc)
		col83 = pf.Column(name = 'FRACEXP', format='E',array = blind2exp)
		col84 = pf.Column(name = 'ERROR', format='E',unit='counts',array = blind2err)
		col85 = pf.Column(name = 'DEADC', format='E',array = blind2de)
		cols8 = pf.ColDefs([col81,col82,col83,col84,col85])
		tbhdu8 = pf.BinTableHDU.from_columns(cols8,name="BLINDDETECTOR3")
		#########################################header##############################################
		head_file=pf.open("%s/me_0_small_g0_0-17.lc"%Net_obspath)
		head_or=head_file[1].header
		hdr0 = pf.Header()
		hdr0['TELESCOP']=head_or['TELESCOP']
		hdr0['INSTRUME']=head_or['INSTRUME']
		hdr0['OBS_MODE']=head_or['OBS_MODE']
		hdr0['OBS_ID']=head_or['OBS_ID']
		hdr0['OBS_MJD']=head_or['MJD-OBS']
		hdr0['OBS_DATA']=head_or['DATE-OBS']
		c_icrs = SkyCoord(ra=head_or['RA_PNT']*u.degree, dec=head_or['DEC_PNT']*u.degree, frame='icrs')
		hdr0['PNT_B']=('%s'%(int(round(c_icrs.galactic.b.degree))),'[deg],Galactic Latitude')
		hdr0['PNT_L']=('%s'%(int(round(c_icrs.galactic.l.degree))),'[deg],Galactic Longitude')
		hdr0['TSTART']=(head_or['TSTART'],'Start Time')
		hdr0['TSTOP']=(head_or['TSTOP'],'Stop Time')
		hdr0['SOFTWARE']=('hxmtsoftV2','Light Curve Tools')
		file=os.stat(r"%s/%s_ub_chi2.dat" % (Midd_obspath, ObsID))
		ftime=time.asctime(time.localtime(file.st_mtime))
		hdr0['TFITTED']=('%s'%ftime,'Fiting Date')
		prhdu = pf.PrimaryHDU(header=hdr0)
		####################################build file##########################################
		tbhdulist2 = pf.HDUList([prhdu,tbhdu1,tbhdu2,tbhdu3,tbhdu4,tbhdu5,tbhdu6,tbhdu7,tbhdu8])
		tbhdulist2.writeto("%s/%s_FITDATA_ME.fits" % (Midd_obspath, ObsID), clobber=True)
		mkdir_try(fProd_obspath)
		os.system('cp %s/%s_FITDATA_ME.fits %s/%s_FITDATA_ME.fits'% (Midd_obspath, ObsID, fProd_obspath, ObsID))
		print "%s: finish"%path


if __name__ == "__main__":
	Wpath='/sharefs/hbkg/user/saina/data294/P0201013/'
	#lclist= os.listdir(Wpath)
	#lclist.sort()
	lclist='P020101300101'
	path=Wpath+lclist+'/ME/'
	ehkfile=Wpath+lclist+"/ehk.fits"
	gtifile=path+'me_gti.fits'
	dataout(Wpath,lclist)
