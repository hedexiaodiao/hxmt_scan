import datetime
import time
import os
import numpy as np
import sys
import smtplib
from multiprocessing import Process
from email.mime.text import MIMEText
import project_rdata
import src_map_auto
#sys.Wpath.append("/sharefs/hbkg/user/saina/tar/P0211007")
#import timing_run
import write_007xml_index
#sys.Wpath.append('/sharefs/hbkg/user/saina/tar/psf_201904')
#import createpha_betax as creatpha


def loadmod(ObsID, Wpath):
	#Wpath='/sharefs/hbkg/user/saina/data294/P0211007/'
	OrgWpath = Wpath + '/Org'  ###
	NetWpath = Wpath + '/Net'  ###
	Program_dir = '/sharefs/hbkg/user/luoqi/HXMT_SCAN/ME/ME_SCAN'###
	Midd_obspath = Wpath + '/Midd/' + ObsID  ###
	if not os.path.exists("%s/%s/me_lc_box2_small.fits"%(NetWpath, ObsID)):
		print ObsID, ': no data'
		return 0
	if os.path.exists( "%s/%s_ub_chi2.dat" %(Midd_obspath,ObsID)):
		print ObsID, ' : already done'
		if not os.path.exists("%s/%s_smap.png" %(Midd_obspath,ObsID)):
			project_rdata.dataout(Wpath, ObsID)
			src_map_auto.src_map(Wpath, ObsID)
		return 0
	else:
		if os.path.exists(Midd_obspath):
			os.system("rm -r %s" % (Midd_obspath))
		os.system("cp %s/test_all_typeSub.py %s/test_all_typeSub.py" % (Program_dir,Midd_obspath))
		if not os.path.exists('%s/config_me_small.xml'%(Midd_obspath)):
			write_007xml_index.write_xml(Wpath, ObsID)
		###os.system("mkdir %s/%s/ME/result/" % (Wpath, ObsID))
		print ObsID, " : Fit start"
		os.system("/hxmt/soft/Develop/anaconda2/bin/python %s/test_all_typeSub.py %s/config_me_small.xml" % (Midd_obspath, Midd_obspath))
		project_rdata.dataout(Wpath, ObsID)
		src_map_auto.src_map(Wpath, ObsID)
		return 1

def doSth(list5,fpath):
	plog = []
	for k in list5[:-5]:
		for i in k:
			p = Process(target = loadmod, args = (i,fpath))
			p.start()
			plog.append(p)
		
		for p in plog:
			p.join()
	
	print "Check END"
	time.sleep(60)



def main(h=0, m=0):
	Program_dir = '/sharefs/hbkg/user/luoqi/HXMT_SCAN/ME/ME_SCAN'
	Wpath = '/sharefs/hbkg/data/SCAN/luoqi/ME'
	NetWpath = Wpath + '/Net'  ###
	while True:
		while True:
			now = datetime.datetime.now()
			#if (now.hour%3)==0 and (now.minute%5)==m:
			if (now.minute%5)==m:
				time.sleep(60)
				break
			time.sleep(60)
		
		stime=time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time()))
		print stime
		#lclist=np.loadtxt("/sharefs/hbkg/user/saina/tar/P0211007/aplist.txt", delimiter='\t', dtype=bytes).astype(str)
		#obin='P0211007'
		if len(sys.argv)<2:
			print "Need the ffn!"
			obin='P0'+str(input("ff ObsID (211007 211009 201013 101294 101295): "))
			print obin
		else:
			print sys.argv[1]
			obin = sys.argv[1]
		'''
		if (obin[-1]=='7')|(obin[-1]=='9')|(obin[-1]=='4')|(obin[-1]=='3'):
			lcpath='/sharefs/hbkg/user/saina/data294/'+obin+'/'
		if obin[-1]=='5':
			lcpath='/sharefs/hbkg/user/saina/data295/'+obin+'/'
		'''
		lcpath = NetWpath###Important, use exist Net files only
		lclist=os.listdir(lcpath)
		lclist.sort()
		os.chdir('%s/lcfit/psf_201904/'%Program_dir)
		os.system('/hxmt/soft/Develop/anaconda2/bin/python creatalpha_betax.py')
		os.chdir(Wpath)
		fitlist=lclist[::-1]###
		fitnum=[]
		Wpath=lcpath
		fob=0
		list5=[[]]*int(1+len(fitlist)/5)###5 ObsID 1 group
		for i in range(len(list5)):
			list5[i]=fitlist[i*5:(i+1)*5]###
			#print list5[i]
		#rlist=[]
		doSth(list5,lcpath)
		for k in list5:
			for ObsID in k:
				Midd_obspath = Wpath + '/Midd/' + ObsID  ###
				if os.path.exists('%s/%s_smap.png'%(Midd_obspath,ObsID))&os.path.exists("%s/%s_ub_chi2.dat"%(Midd_obspath,ObsID)):
					fitnum=np.append(fitnum,ObsID)
					print ObsID,' : Fit end'
				else:
					print ObsID,' : Fit error!'
		
		'''for ObsID in fitlist:
			fob=loadmod(ObsID,Wpath)
			if os.Wpath.exists('%s/%s/ME/result/%s_smap.png'%(Wpath,ObsID,ObsID))&(fob>0):
				fitnum=np.append(fitnum,ObsID)
				print ObsID,' : Fit end'
			fob=0'''
		
		titleld = '2--New fit list: %s'%stime
		if len(fitnum)>0:
			contentld = '%s'%fitnum
		
		print 'fit end'



main(0,0)
