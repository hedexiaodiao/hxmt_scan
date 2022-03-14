#!/hxmt/soft/Develop/anaconda2/bin/python
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.io import fits as pf
import pandas as pd
from astropy.time import Time
from scipy.optimize import curve_fit
import sys
import os
import glob
from decimal import *
def smooth(x,window_len=11,window='hanning'):
	if x.ndim != 1:
		raise ValueError, "smooth only accepts 1 dimension arrays."
	
	if x.size < window_len:
		raise ValueError, "Input vector needs to be bigger than window size."
	
	if window_len<3:
		return x
	
	x_right = x[-(window_len-1+50):] - (np.sum(x[-60:-45])/15 - np.sum(x[-15:])/15)
	x_left = x[:window_len-1+50] - (np.sum(x[45:60])/15 - np.sum(x[0:15])/15)
	s=np.r_[x_left,x,x_right]
	if window == 'flat': #moving average
		w=np.ones(window_len,'d')
	else:
		w=eval('np.'+window+'(window_len)')
	smtcut = (window_len-1)/2 + 50
	y=np.convolve(w/w.sum(),s,mode='valid')[smtcut:-smtcut]
	#print "smooth::: len smooth return:",(y.shape)[0]
	return y

def smooth_bkg(x,window_len=11,window='hanning'):
	x_right = x[-1]*2 - x[-1:-(window_len-1+50):-1]
	x_left = x[0]*2 - x[(window_len-1+50):1:-1]
	s=np.r_[x_left,x,x_right]
	if window == 'flat':
		w=np.ones(window_len,'d')
	else:
		w=eval('np.'+window+'(window_len)')
	smtcut = (window_len-1)/2 + 49
	y=np.convolve(w/w.sum(),s,mode='valid')[smtcut:-smtcut]
	return y


'''
def f_fit(x,r,R,a,b,c,A,B,C,m,n):
	return 0.01*m * (1*r*np.sin(0.0004*a*x+b)+c)*(1*R*np.sin(0.0020*A*x+B)+C) + 10*n#+e*np.cos(0.001*f*x)+g
def f_show(x,p_fit):  
	r,R,a,b,c,A,B,C,m,n=p_fit.tolist()
	#print r,R,a,b,c,A,B,C,m,n
	#x = x + 0.001*e*np.cos(0.001*f*x)
	return 0.01*m * (1*r*np.sin(0.0004*a*x+b)+c)*(1*R*np.sin(0.0020*A*x+B)+C) + 10*n#+e*np.cos(0.001*f*x)+g
'''

def poly_fit(x,y):
	polmin = np.poly1d(np.polyfit(x, y, 2))
	while 1:
		fragres = polmin(x) - y
		if (fragres < -3).sum()>0:
			y = y[fragres!=fragres.min()]
			x = x[fragres!=fragres.min()]
			polmin = np.poly1d(np.polyfit(x, y, 2))
		else:
			break
	while 1:
		fragres = y - polmin(x)
		if (fragres < -5).sum()>0:
			y = y[fragres!=fragres.min()]
			x = x[fragres!=fragres.min()]
			polmin = np.poly1d(np.polyfit(x, y, 2))
		else:
			break

def get_bkg(lcfile):
	hdu = pf.open(lcfile)
	tb = hdu[1].data
	x = tb.field(0)
	y = tb.field(1)
	err = tb.field(3)
	xr = x#[y>0]
	yr = y#[y>0]
	err = err#[y>0]
	Ygti = ((xr[1:]-xr[:-1])!=1)
	c=np.nonzero(Ygti)[0]
	pot=np.r_[0,c+1,Ygti.shape[0]]
	#print pot,"pot::",i
	par = []
	for i in range(len(pot)-1): 
		print pot,"pot::",i
		try:
			yr1 = yr[pot[i]:pot[i+1]][2:-2]
			xr1 = xr[pot[i]:pot[i+1]][2:-2]
			err1 = err[pot[i]:pot[i+1]][2:-2]
		except IndexError:
			raw_input("Indexerror ,recompose you code!")
			yr1 = yr[pot[i]:pot[i+1]]
			xr1 = xr[pot[i]:pot[i+1]]
			err1 = err[pot[i]:pot[i+1]]
		print ':::::',len(xr1),len(yr1)
		plt.plot(xr1,yr1,'r-',alpha=0.5)
		if yr1.shape[0]>=600:
			ysmooth = smooth(yr1,window_len=11)
			yr1[(yr1-ysmooth)>60] = ysmooth[(yr1-ysmooth)>60]
			ysmooth = smooth(yr1,window_len=5,window='blackman')
			fraglen = 30
			startnum = 300/fraglen
			nmb = yr1.shape[0]/fraglen+2*startnum
			steplen = 15
			lcdic = {}
			for idx in range(nmb):
				stat = (idx-startnum) * fraglen  if (idx-startnum)*fraglen > 0 else 0
				stop = (idx-startnum) * fraglen + 600 if (idx-startnum) * fraglen + 600 < yr1.shape[0] else yr1.shape[0]
				#print stat,stop,yr1.shape,idx
				if stat>=yr1.shape[0]-300:
					stat=yr1.shape[0]-300
					break
				frag = ysmooth[stat:stop]
				fragmin = ysmooth[stat:stop]
				#print fragmin.shape
				tmmin = xr1[stat:stop]
				polmin = np.poly1d(np.polyfit(tmmin-xr1[0], fragmin, 2))
				while 1:
					fragres = polmin(tmmin-xr1[0]) - fragmin
					if (fragres < -10).sum()>0:
						fragmin = fragmin[fragres!=fragres.min()]
						tmmin = tmmin[fragres!=fragres.min()]
						polmin = np.poly1d(np.polyfit(tmmin-xr1[0], fragmin, 2))
					else:
						break
				while 1:
					fragres = fragmin - polmin(tmmin-xr1[0])
					if (fragres < -10).sum()>0:
						fragmin = fragmin[fragres!=fragres.min()]
						tmmin = tmmin[fragres!=fragres.min()]
						polmin = np.poly1d(np.polyfit(tmmin-xr1[0], fragmin, 2))
					else:
						break
				while 1:
					fragres = polmin(tmmin-xr1[0]) - fragmin
					if (fragres < -5).sum()>0:
						fragmin = fragmin[fragres!=fragres.min()]
						tmmin = tmmin[fragres!=fragres.min()]
						polmin = np.poly1d(np.polyfit(tmmin-xr1[0], fragmin, 2))
					else:
						break
				while 1:
					fragres = fragmin - polmin(tmmin-xr1[0])
					if (fragres < -5).sum()>0:
						fragmin = fragmin[fragres!=fragres.min()]
						tmmin = tmmin[fragres!=fragres.min()]
						polmin = np.poly1d(np.polyfit(tmmin-xr1[0], fragmin, 2))
					else:
						break
				while 1:
					fragres = polmin(tmmin-xr1[0]) - fragmin
					if (fragres < -3).sum()>0:
						fragmin = fragmin[fragres!=fragres.min()]
						tmmin = tmmin[fragres!=fragres.min()]
						polmin = np.poly1d(np.polyfit(tmmin-xr1[0], fragmin, 2))
					else:
						break
				
				for ix in range(len(xr1[stat:stop])):
					if lcdic.has_key('%.8f'%(xr1[stat:stop][ix])):
						lcdic['%.8f'%(xr1[stat:stop][ix])].append(polmin(xr1[stat:stop][ix]-xr1[0]))
					else:
						lcdic['%.8f'%(xr1[stat:stop][ix])]= [polmin(xr1[stat:stop][ix]-xr1[0])]
				plt.plot(xr1[stat:stop],polmin(xr1[stat:stop]-xr1[0]),'g-',lw=1)
				#plt.plot(xr1[stat:stop],polmax(xr1[stat:stop]-xr1[0]),'g-',lw=1)
				plt.plot(tmmin,fragmin,'b.',ms=8)
				#plt.plot(tmmax,fragmax,'r.',ms=8)
			for key,vals in lcdic.items():
				lcdic[key]=np.median(vals)
			tptm = sorted(lcdic.keys())
			tplc = [lcdic[tmid] for tmid in tptm]
			tptms = [np.float64(tmid) for tmid in tptm]
			ystp = smooth_bkg(np.array(tplc),window_len=151,window='hamming')
			#ystp = smooth_bkg(np.array(ystp),window_len=151)
			#ystp = smooth_bkg(np.array(ystp),window_len=151)
			for step in range(0+steplen,len(xr1),steplen):
				frag = ysmooth[step:step+steplen]
				tp = frag.min()
				timeloc = xr1[step:step+steplen][frag==tp][0]
				#print tp,timeloc,tptms,len(tptms),len(ystp),len(xr1)
				print (np.array(tptms)==timeloc).sum()
				if np.array(ystp)[np.array(tptms, dtype='i')==int(timeloc)][0]>tp:
					headcut = timeloc - xr1[0]
					break
			for step in range(len(xr1)-1-steplen,0,-steplen):
				frag = ysmooth[step - steplen:step]
				tp = frag.min()
				timeloc = xr1[step - steplen:step][frag==tp][-1]
				#print tp,timeloc,len(tptms),len(ystp),len(xr1)
				if np.array(ystp)[np.array(tptms, dtype='i')==int(timeloc)][0]>tp:
					endcut = timeloc - xr1[-1]
					break
			print headcut,endcut
			plt.plot(tptms,tplc,'b-',lw=1.5)
			plt.plot(tptms,ystp,'y-',lw=1.5)
			#plt.show()
			bkg1 = ystp[int(headcut+60):int(endcut-60)]
			xr1 = tptms[int(headcut+60):int(endcut-60)]
			yr1 = yr1[int(headcut+60):int(endcut-60)]
			err1 = err1[int(headcut+60):int(endcut-60)]
			if i ==0:
				bkg = bkg1
				xr_sel = xr1
				yr_sel = yr1
				err_sel = err1
			if i >0:
				try:
					bkg = np.r_[bkg,bkg1]
					xr_sel = np.r_[xr_sel,xr1]
					yr_sel = np.r_[yr_sel,yr1]
					err_sel = np.r_[err_sel,err1]
				except UnboundLocalError:
					bkg = bkg1
					xr_sel = xr1
					yr_sel = yr1
					err_sel = err1
		#plt.close('all')
	mask = np.in1d(x,xr_sel)
	yr_sel = y[mask]
	hdu.close()
	del hdu
	return xr_sel,yr_sel,bkg,err_sel

def poly_bkg(NetWpath, Wpath, ObsID):
	#key = sys.argv[1]
	#inpath = "/sharefs/hbkg/user/saina/data295/P0101295/" + key + "/"
	#outpath = "/sharefs/hbkg/user/saina/data295/P0101295/"
	inpath = NetWpath + '/' + ObsID ###
	outpath = inpath
	slcfile_list=glob.glob(r'%s/*small.lc'%(inpath))
	#blcfile_list=glob.glob(r'%s/*blind_g0*.lc'%(inpath))
	slcfile_list.sort()
	#blcfile_list.sort()
	xr0_small,yr0_small,bkg0_small,error0 = get_bkg(slcfile_list[0])
	xr0 = np.r_[xr0_small]
	yr0 = np.r_[yr0_small]
	bkg0 = np.r_[bkg0_small]
	col1 = pf.Column(name = 'Time', format='D',array = xr0)
	col2 = pf.Column(name = 'Counts', format='D',array = yr0-bkg0)
	col3 = pf.Column(name = 'Stat_err', format='D',array = error0)
	col4 = pf.Column(name = 'Bkg', format='D',array = bkg0)
	cols = pf.ColDefs([col1,col2,col3])
	tbhdu = pf.BinTableHDU.from_columns(cols)
	tbhdu.writeto(outpath+"/me_lc_box0_small_bkg.fits",clobber=True)
	cols2 = pf.ColDefs([col1,col4])
	tbhdu2 = pf.BinTableHDU.from_columns(cols2)
	tbhdu2.writeto(outpath+"/bkg0.fits",clobber=True)
	print "me 0 done"
	
	xr1_small,yr1_small,bkg1_small,error1 = get_bkg(slcfile_list[1])
	xr1 = np.r_[xr1_small]
	yr1 = np.r_[yr1_small]
	bkg1 = np.r_[bkg1_small]
	col1 = pf.Column(name = 'Time', format='D',array = xr1)
	col2 = pf.Column(name = 'Counts', format='D',array = yr1-bkg1)
	col3 = pf.Column(name = 'Stat_err', format='D',array = error1)
	col4 = pf.Column(name = 'Bkg', format='D',array = bkg1)
	cols = pf.ColDefs([col1,col2,col3])
	tbhdu = pf.BinTableHDU.from_columns(cols)
	tbhdu.writeto(outpath+"/me_lc_box1_small_bkg.fits",clobber=True)
	cols2 = pf.ColDefs([col1,col4])
	tbhdu2 = pf.BinTableHDU.from_columns(cols2)
	tbhdu2.writeto(outpath+"/bkg1.fits",clobber=True)
	print "me 1 done"
	
	xr2_small,yr2_small,bkg2_small,error2 = get_bkg(slcfile_list[2])
	xr2 = np.r_[xr2_small]
	yr2 = np.r_[yr2_small]
	bkg2 = np.r_[bkg2_small]
	col1 = pf.Column(name = 'Time', format='D',array = xr2)
	col2 = pf.Column(name = 'Counts', format='D',array = yr2-bkg2)
	col3 = pf.Column(name = 'Stat_err', format='D',array = error2)
	col4 = pf.Column(name = 'Bkg', format='D',array = bkg2)
	cols = pf.ColDefs([col1,col2,col3])
	tbhdu = pf.BinTableHDU.from_columns(cols)
	tbhdu.writeto(outpath+"/me_lc_box2_small_bkg.fits",clobber=True)
	cols2 = pf.ColDefs([col1,col4])
	tbhdu2 = pf.BinTableHDU.from_columns(cols2)
	tbhdu2.writeto(outpath+"/bkg2.fits",clobber=True)
	print "me 2 done"
	
	fig,axe = plt.subplots(3,sharex=True, sharey=False,figsize=(22,16))
	fig.subplots_adjust(hspace = 0.05,top=0.9,bottom=0.1,left=0.1,right=0.95)
	l0=axe[0].plot(xr0,yr0,c='b',lw=0.5,label=r'$60^o $')
	l1=axe[1].plot(xr1,yr1,c='b',lw=0.5,label=r'$0^o$')
	l2=axe[2].plot(xr2,yr2,c='b',lw=0.5,label=r'$-60^o$')
	l3=axe[0].plot(xr0,bkg0,c='r',lw=0.5)
	l4=axe[1].plot(xr1,bkg1,c='r',lw=0.5)
	l5=axe[2].plot(xr2,bkg2,c='r',lw=0.5)
	l6=axe[0].plot(xr0,yr0-bkg0,c='k',lw=0.5)
	l7=axe[1].plot(xr1,yr1-bkg1,c='k',lw=0.5)
	l8=axe[2].plot(xr2,yr2-bkg2,c='k',lw=0.5)
	axe[1].set_ylabel('counts/sec',fontsize=12)
	
	for i in range(3):
		axe[i].grid()
		axe[i].legend(loc='upper right')
		axe[i].legend(labels=['data','blackground','residual'], loc='upper left')
	
	#fig.legend((l0,l3,l6),('diata','blackground','residual'),'center right')
	
	#axe[0].set_title(key[:13]+"_he")
	plt.xlabel('Time /s',fontsize = 12)
	plt.savefig(outpath+"/me_lc_small_bkg.png",dpi=400)
	#plt.show()
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
	path = '/hxmt/work/HXMT-DATA/1L/A%s/%s/' % (ObsID[1:3], ObsID[:-5])

	###Wpath='/sharefs/hbkg/user/saina/data294/P0201013/'
	#lclist= os.listdir(Wpath)
	#lclist.sort()
	###lclist='P020101300201'
	#Wpath=Wpath+lclist
	#ehkfile=Wpath+lclist+"/ehk.fits"
	#gtifile=Wpath+'me_gti.fits'
	poly_bkg(Wpath,Wpath,ObsID)
