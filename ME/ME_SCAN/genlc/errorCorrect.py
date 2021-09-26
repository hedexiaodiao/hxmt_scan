#!/hxmt/soft/Develop/anaconda2/bin/python
import os
import os.path
import numpy as np
import astropy.io.fits as pf
import glob
import sys

def errorCorrect(Wpath, ObsID):
	OrgWpath = Wpath + '/Org'  ###
	Org_obspath = "%s/%s" % (OrgWpath, ObsID)  ###
	Org_screenpath = "%s/%s" % (OrgWpath, 'Screen')  ###
	NetWpath = Wpath + '/Net'  ###
	Wscreenfile = "%s/%s_me_screen.fits" % (Org_screenpath, ObsID)  ###
	Wdeadfile = "%s/%s_me_dead.fits" % (Org_obspath, ObsID)  ###dead time file
	'''Ppath='/sharefs/hbkg/user/saina/pfiles_small/'
	#if os.Wpath.exists(Ppath+"/pfiles_"+ObsID):
		#os.system("rm -r %s/pfiles_%s"%(Ppath,ObsID))
	os.system("mkdir %s/pfiles_%s"%(Ppath,ObsID))
	os.environ["PFILES"]="%s/pfiles_%s;/home/hxmt/hxmtsoft2/soft/install/x86_64-unknown-linux-gnu-libc2.12/syspfiles"%(Ppath,ObsID)
	#os.system("echo $PFILES")'''
	me_lc_cmd = '''for box in 0 1 2; do
			case $box in
			"0")
			small="0-7,11-17"
			big="8-9"
			blind="10"
			;;
			"1")
			small="18-25,29-35"
			big="26-27"
			blind="28"
			;;
			"2")
			small="36-43,47-53"
			big="44-45"
			blind="46"
			;;
			esac
	melcgen evtfile=%s deadfile=%s outfile=%s/%s/%s_me_${box}_smallErr userdetid="${small}" starttime=0 stoptime=0 minPI=68 maxPI=631 binsize=1 deadcorr=no
	melcgen evtfile=%s deadfile=%s outfile=%s/%s/%s_me_${box}_blindErr userdetid="${blind}" starttime=0 stoptime=0 minPI=68 maxPI=631 binsize=1 deadcorr=no
		done'''%(Wscreenfile, Wdeadfile, NetWpath, ObsID, ObsID, Wscreenfile, Wdeadfile, NetWpath, ObsID, ObsID)
	'''os.system(me_lc_cmd)
	os.system('rm -r %s/pfiles_%s'%(Ppath,ObsID))'''
	slcfile_list=glob.glob(r'%s/%s/*small_g0*.lc' % (NetWpath, ObsID))
	blcfile_list=glob.glob(r'%s/%s/*blind_g0*.lc' % (NetWpath, ObsID))
	#selc_list=glob.glob(r'%s/%s/*/*smallErr_g0*.lc'%(Wpath,ObsID))
	#belc_list=glob.glob(r'%s/%s/*/*blindErr_g0*.lc'%(Wpath,ObsID))
	slcfile_list.sort()
	blcfile_list.sort()
	#selc_list.sort()
	#belc_list.sort()
	#namelist=['0-17','18-35','36-53']
	#namelist_b=['10','28','46']
	for num in range(3):
		data1=pf.open(slcfile_list[num])
		data_or=data1[1].data
		#data2=pf.open(selc_list[ObsID])
		#data_de=data2[1].data
		#mask1=np.in1d(data_or.field(0),data_de.field(0))
		#mask2=np.in1d(data_de.field(0),data_or.field(0))
		count1=data_or.field(1)
		#[mask1]
		count2=count1
		#count2=data_de.field(1)[mask2]
		#error_de=data_de.field(3)
		error_cal=[0]*len(count2)
		for i in range(len(count2)):
			if count2[i]!=0:
				error_cal[i]=(1+np.sqrt(count2[i]+0.75))
				#*(count1[i]/count2[i])
			else:
				error_cal[i]=1+np.sqrt(count2[i]+0.75)
		#print "step1"
		col1 = pf.Column(name = 'TIME', format='1D',unit='s',array = data_or.field(0))#[mask1])
		col2 = pf.Column(name = 'COUNTS', format='1E',unit='counts',array = data_or.field(1))#[mask1])
		col3 = pf.Column(name = 'FRACEXP', format='1E',array = data_or.field(2))#[mask1])
		col4 = pf.Column(name = 'ERROR', format='1E',unit='counts',array = error_cal)
		lccols = pf.ColDefs([col1,col2,col3,col4])
		lctbhdu = pf.BinTableHDU.from_columns(lccols,name='COUNTS')
		lenfile = len(data1)
		data1hdu =data1[lenfile-1]
		tbhdu=[[] for j in range(lenfile-3)]
		for i in range(lenfile-3):
			tbhdu[i] = data1[i+2]
		prihdu = data1[0]
		list1=[prihdu,lctbhdu]
		for k in range(len(tbhdu)):
			list1=np.append(list1,tbhdu[k])
		list1=np.append(list1,data1hdu)
		list1=list(list1)
		tbhdulist = pf.HDUList(list1)
		tbhdulist.writeto("%s/%s/me_%s_small.lc" % (NetWpath, ObsID, num), clobber=True)
		##################################blind############################################
		datab1=pf.open(blcfile_list[num])
		data_orb=datab1[1].data
		#datab2=pf.open(belc_list[ObsID])
		#data_deb=datab2[1].data
		#maskb1=np.in1d(data_orb.field(0),data_deb.field(0))
		#maskb2=np.in1d(data_deb.field(0),data_orb.field(0))
		countb1=data_orb.field(1)
		#[maskb1]
		countb2=countb1
		#countb2=data_deb.field(1)[maskb2]
		#error_de=data_de.field(3)
		error_calb=[0]*len(countb2)
		for j in range(len(countb2)):
			if countb2[j]!=0:
				error_calb[j]=(1+np.sqrt(countb2[j]+0.75))
				#*(countb1[j]/countb2[j])
			else:
				error_calb[j]=1+np.sqrt(countb2[j]+0.75)
		#print "step1"
		colb1 = pf.Column(name = 'TIME', format='1D',unit='s',array = data_orb.field(0))#[maskb1])
		colb2 = pf.Column(name = 'COUNTS', format='1E',unit='counts',array = data_orb.field(1))#[maskb1])
		colb3 = pf.Column(name = 'FRACEXP', format='1E',array = data_orb.field(2))#[maskb1])
		colb4 = pf.Column(name = 'ERROR', format='1E',unit='counts',array = error_calb)
		colb5 = pf.Column(name = 'DEADC', format='1E',array = data_orb.field(4))#[maskb1])
		lccolbs = pf.ColDefs([colb1,colb2,colb3,colb4,colb5])
		lctbbhdu = pf.BinTableHDU.from_columns(lccolbs,name='COUNTS')
		lenfileb = len(datab1)
		datab1hdu =datab1[lenfileb-1]
		tbbhdu=[[] for j in range(lenfileb-3)]
		for i in range(lenfileb-3):
			tbbhdu[i] = datab1[i+2]
		pribhdu = datab1[0]
		listb1=[pribhdu,lctbbhdu]
		for k in range(len(tbbhdu)):
			listb1=np.append(listb1,tbbhdu[k])
		listb1=np.append(listb1,datab1hdu)
		listb1=list(listb1)
		tbhdublist = pf.HDUList(listb1)
		tbhdublist.writeto("%s/%s/me_%s_blind.lc" % (NetWpath, ObsID, num), clobber=True)
'''
	data1=pf.open("%s/%s/ME/me_1_small_g0_18-35.lc"%(Wpath,ObsID))
	data_or=data1[1].data
	data2=pf.open("%s/%s/ME/me_1_smallEr_g0_18-35.lc"%(Wpath,ObsID))
	data_de=data2[1].data
	error_or=data_or.field(3)
	error_cal=1+np.sqrt(error_or**2+0.75)
	col1 = pf.Column(name = 'TIME', format='1D',unit='s',array = data_de.field(0))
	col2 = pf.Column(name = 'COUNTS', format='1E',unit='counts',array = data_de.field(1))
	col3 = pf.Column(name = 'FRACEXP', format='1E',array = data_de.field(2))
	col4 = pf.Column(name = 'ERROR', format='1E',unit='counts',array = error_cal)
	lccols = pf.ColDefs([col1,col2,col3,col4])
	lctbhdu = pf.BinTableHDU.from_columns(lccols,name='COUNTS')
	lenfile = len(data2)
	data2hdu =data2[lenfile-1]
	tbhdu=[[] for j in range(lenfile-3)]
	for i in range(lenfile-3):
		tbhdu[i] = data2[i+2]
	prihdu = data2[0]
	list1=[prihdu,lctbhdu]
	for k in range(len(tbhdu)):
		list1=np.append(list1,tbhdu[k])
	list1=np.append(list1,data2hdu)
	list1=list(list1)
	tbhdulist = pf.HDUList(list1)
	tbhdulist.writeto("me_1_small.lc",clobber=True)

	data1=pf.open("%s/%s/ME/me_2_small_g0_36-53.lc"%(Wpath,ObsID))
	data_or=data1[1].data
	data2=pf.open("%s/%s/ME/me_2_smallErr_g0_36-53.lc"%(Wpath,ObsID))
	data_de=data2[1].data
	error_or=data_or.field(3)
	error_cal=1+np.sqrt(error_or**2+0.75)
	col1 = pf.Column(name = 'TIME', format='1D',unit='s',array = data_de.field(0))
	col2 = pf.Column(name = 'COUNTS', format='1E',unit='counts',array = data_de.field(1))
	col3 = pf.Column(name = 'FRACEXP', format='1E',array = data_de.field(2))
	col4 = pf.Column(name = 'ERROR', format='1E',unit='counts',array = error_cal)
	lccols = pf.ColDefs([col1,col2,col3,col4])
	lctbhdu = pf.BinTableHDU.from_columns(lccols,name='COUNTS')
	lenfile = len(data2)
	data2hdu =data2[lenfile-1]
	tbhdu=[[] for j in range(lenfile-3)]
	for i in range(lenfile-3):
		tbhdu[i] = data2[i+2]
	prihdu = data2[0]
	list1=[prihdu,lctbhdu]
	for k in range(len(tbhdu)):
		list1=np.append(list1,tbhdu[k])
	list1=np.append(list1,data2hdu)
	list1=list(list1)
	tbhdulist = pf.HDUList(list1)
	tbhdulist.writeto("me_2_small.lc",clobber=True)
'''
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

	###Wpath='/sharefs/hbkg/user/saina/data294/P0211007/'
	#lclist= os.listdir(Wpath)
	#lclist.sort()
	###lclist='P021100703201'
	#Wpath=Wpath+lclist+'/ME/'
	errorCorrect(Wpath,ObsID)