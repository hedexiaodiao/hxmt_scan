#!/hxmt/home/lick/soft/anaconda2/bin/python
'''
Copyright 2016, HXMT Ground System
All right reserved
File Name: loadmod.py for sample

.......

Author: lick
Version: 0.1
Date:  2016-12-11
Email: lick@ihep.ac.cn
......
'''
from astropy.io import fits as pf
import sys
import os
import time
import numpy as np
import math
def cpha(crf,fake):
	fname = crf.split('/')[-2][9:13]
	temppath = fake[:-8]
	print('cpha para, crf, fake, fname:', crf, fake, fname)
	print('temppath:',temppath)
	hd = pf.open(crf)
	tb = hd[1].data
	t = tb.field(0)
	t0  = t[0]
	cr2 = tb.field(1)
	err2 = tb.field(2)
	os.system('pwd')
	ran = str(time.time())
	f=open("%s%stemppha%s.txt"%(temppath,fname,ran),"w")
	f2=open("%s%stempbad%s.dat"%(temppath,fname,ran),"w")
	for i,r in enumerate(cr2):
		tst = i 
		tsp = i+1
		r = r
		if i ==0:
			if r <1:
				tl = str(i)+" "
				f2.writelines(tl)
		if i >0 and i <len(cr2)-1:
			if r<1:
				if cr2[i-1]>1:
					tl2 = str(i)+" "
					f2.writelines(tl2)
			if r>1:
				if cr2[i-1]<1:
					tl3 = str(i-1)+"\n"
					f2.writelines(tl3)
		if i == len(cr2)-1:
			if r<1:
				if cr2[i-1]<1:
					tl4 =str(i)+"\n"
					f2.writelines(tl4)
				if cr2[i-1]>=1:
					tl4 =str(i)+" "+str(i)+"\n"
					f2.writelines(tl4)
		re = err2[i]
		line = str(tst)+" "+str(tsp)+" "+str(r)+" "+str(re)+"\n"
		f.writelines(line)
	for i,r in enumerate(err2):
		tst = i 
		tsp = i+1
		r = r
		if i ==0:
			if r <0.1:
				tl = str(i)+" "
				f2.writelines(tl)
		if i >0 and i <len(err2)-1:
			if r<0.1:
				if err2[i-1]>0.1:
					tl2 = str(i)+" "
					f2.writelines(tl2)
			if r>0.1:
				if err2[i-1]<0.1:
					tl3 = str(i-1)+"\n"
					f2.writelines(tl3)
		if i == len(err2)-1:
			if r<0.1:
				if err2[i-1]<0.1:
					tl4 =str(i)+"\n"
					f2.writelines(tl4)
				if err2[i-1]>=0.1:
					tl4 =str(i)+" "+str(i)+"\n"
					f2.writelines(tl4)
		#print line
	f.close()
	f2.close()
	os.system('ls %s%stemppha%s.txt'%(temppath,fname,ran))
	#os.environ[
	#	"PFILES"] = "%s/pfiles_%s_load;/hxmt/soft/Astro/heasoft-6.27.2/x86_64-pc-linux-gnu-libc2.17/syspfiles" % (
	#Ppath, obid)
	###commond = 'export PFILES="%s/pfiles_%s_load;/hxmt/soft/Astro/heasoft-6.27.2/x86_64-pc-linux-gnu-libc2.17/syspfiles";export HEADASNOQUERY=;export HEADASPROMPT=/dev/null;flx2xsp %s%stemppha%s.txt %s.pha %s.rsp' % (
	###Ppath, 'P010129500103', temppath, fname, ran, fake, fake)
	commond = 'export HEADASNOQUERY=;export HEADASPROMPT=/dev/null;flx2xsp %s%stemppha%s.txt %s.pha %s.rsp' % (temppath, fname, ran, fake, fake)
	print(commond)
	os.system(commond)###
	os.system('rm %s%stemppha%s.txt'%(temppath,fname,ran))###
	commond = 'grppha %s.pha !%s.grp "bad %s%stempbad%s.dat&exit"'%(fake,fake,temppath,fname,ran)###
	print(commond)
	os.system(commond)
	#f3=open('%s%stempbad%s.dat'%(crf[:-51],fname,ran))
	#print f3.readlines()
	os.system('rm %s%stempbad%s.dat'%(temppath,fname,ran))
	print("generation of the faked pha file and rsp file : success")

if __name__ == "__main__":
	#print "test"
	Program_dir = '/sharefs/hbkg/user/luoqi/HXMT_SCAN/ME/ME_SCAN'
	Wpath = "/sharefs/hbkg/data/SCAN/luoqi/ME"  ###'/sharefs/hbkg/user/saina/data294/P0201013/'
	MiddWpath = Wpath + '/Midd'  ###
	NetWpath = Wpath + '/Net'  ###

	lclist= os.listdir(NetWpath)
	lclist.sort()
	#lclist=['P010129411301']
	Ppath='%s/pfiles_small'%Program_dir###'/sharefs/hbkg/user/saina/pfiles_small/'
	for obid in lclist:
        #        if os.path.exists('%s/%s/me_3box2.grp'%(MiddWpath,obid)):
        #                continue
		#if not os.path.exists('%s/%s/me_lc_box2_small_cut.fits'%(NetWpath,obid)):
		#	continue
		os.system("mkdir %s/pfiles_%s_load"%(Ppath,obid))
		###os.environ["PFILES"]="%s/pfiles_%s_load;/home/hxmt/hxmtsoft2/soft/install/x86_64-unknown-linux-gnu-libc2.12/syspfiles"%(Ppath,obid)
		os.environ[
			"PFILES"] = "%s/pfiles_%s_load;/hxmt/soft/Astro/heasoft-6.27.2/x86_64-pc-linux-gnu-libc2.17/syspfiles" % (
		Ppath, obid)
		for i in range(3):
			patha=NetWpath+'/'+obid+'/me_lc_box%s_small_cut.fits'%i
			pathb=MiddWpath+'/'+obid+'/me_3box%s'%i
			###patha='../../Net/%s/me_lc_box%s_small_cut.fits'%(obid,i)
			###pathb='./me_3box%s'%i
			os.chdir(MiddWpath+'/'+obid)
			cpha(patha,pathb)
			#rpatha=Wpath+ObsID+'/ME/result/me_lc_box%s_small_res.fits'%i
			#rpathb=Wpath+ObsID+'/me3boxr%s'%i
			#cpha(rpatha,rpathb)
		os.system("rm -r %s/pfiles_%s_load"%(Ppath,obid))
		print(obid,': finish')
