#!/hxmt/soft/Develop/anaconda2/bin/python

import os
import glob
import os.path
import re
import sys
import matplotlib
matplotlib.use('Agg')
import astropy.io.fits as pf
import numpy as np
import regti_v2
import errorCorrect
import gen_lc_me
import mescreen_newbd
import nbkghesmtpoly2
import fix_eci_1N
import check_reatt_auto
os.system("source /hxmt/home/hxmtsoft2/hxmtsoft_v2.sh")
def merun_v2(path,Wpath,num):
	obid = str(num[:-2])
	#Wpath = "/hxmt/work/HXMT-DATA/1L/A02/P0211007/"+ObsID+"/"
	Npath = path[:22]+'N/'+path[28:]
	#Wpath = "/sharefs/hbkg/user/saina/data294/P0211007/"
	print path,Wpath
	Ppath='/sharefs/hbkg/user/saina/pfiles_small/'
	if os.path.exists(Ppath+"/pfiles_"+num):
		os.system("rm -r %s/pfiles_%s"%(Ppath,num))
	os.system("mkdir %s/pfiles_%s"%(Ppath,num))
	os.environ["PFILES"]="%s/pfiles_%s;/home/hxmt/hxmtsoft2/soft/install/x86_64-unknown-linux-gnu-libc2.12/syspfiles"%(Ppath,num)
	os.system("echo $PFILES")
	######################get all the files and folders##########################################################
	attfile_list = glob.iglob(r'%s/%s/*/*_Att_*'%(path,obid))
	meevtfile_list = glob.glob(r'%s/%s/%s*/*/*ME-Evt*'%(path,obid,num))
	metempfile_list = glob.glob(r'%s/%s/%s*/*/*ME-TH*'%(path,obid,num))
	orbitfile_list = glob.glob(r'%s/%s/*/*_Orbit_*'%(path,obid))
	ehkfile_list = glob.iglob(r'%s/%s/*/*EHK*'%(path,obid))
	nattfile_list = glob.iglob(r'%s/%s/*/*_Att_*'%(Npath,obid))
	norbitfile_list = glob.glob(r'%s/%s/*/*_Orbit_*'%(Npath,obid))
	print "%s: files loaded"%obid
	flag1=nflag1=0
	flag2=nflag2=0
	flag3=nflag3=0
	flag4=nflag4=0
	flag5=nflag5=0
	for i in attfile_list:
		if re.findall(r"VU",i)==[]:
			attfile = i
			flag1=1
	for i in meevtfile_list:
		if re.findall(r"VU",i)==[]:
			meevtfile = i
			flag2=1
	for i in metempfile_list:
		if re.findall(r"VU",i)==[]:
			metempfile = i
			flag3=1
	for i in orbitfile_list:
		if re.findall(r"VU",i)==[]:
			orbitfile = i
			flag4=1
	for i in ehkfile_list:
		if re.findall(r"VU",i)==[]:
			ehkfile = i
			flag5=1
	for i in nattfile_list:
		if re.findall(r"VU",i)==[]:
			nattfile = i
			nflag1=1
	for i in norbitfile_list:
		if re.findall(r"VU",i)==[]:
			norbitfile = i
			nflag4=1
	if (flag1!=1)|(flag2!=1)|(flag3!=1)|(flag4!=1):
		print "%s %s: 1st -- Files error"%(obid,num)
	if not os.path.exists("%s/%s/ME"%(Wpath,num)):
		os.makedirs("%s/%s/ME"%(Wpath,num))
	try:
		print orbitfile
	except NameError:
		print "%s %s can not be calculated: orbit"%(obid,num)
		if (nflag4==1)&(nflag1==1):
			print "%s %s 1N orbit test"%(obid,num)
			attfile=nattfile
			fix_eci_1N.fix_eci(norbitfile,"%s/%s/%s_Orbit_1N.fits"%(Wpath,num,num))
			orbitfile="%s/%s/%s_Orbit_1N.fits"%(Wpath,num,num)
		else:
			return 0
	try:
		print attfile
	except NameError:
		print "%s can not be calculated: att"%num
		return 0
	try:
		print meevtfile
	except NameError:
		print "%s can not be calculated: evt"%num
		return 0
	try:
		print metempfile
	except NameError:
		print "%s can not be calculated: th"%num
		return 0
	try:
		print ehkfile
	except NameError:
		print "No ehk file, need calculate by soft"
		flag5=0
	print "---------------------------------------------------------------------------------------"
	###################copy files to the floders##############################################################
	os.system("cp %s %s/%s/att.fits"%(attfile,Wpath,num))
	#os.system("rm %s/%s/ME/gtiv2.fits"%(Wpath,num))
	evt_im = pf.open("%s"%(meevtfile))
	ename = evt_im
	ra_centre = ename[1].header['RA_PNT']
	dec_centre = ename[1].header['DEC_PNT']
	start_time = ename[1].header['TSTART']
	stop_time = ename[1].header['TSTOP']
	print "centre: ",ra_centre,dec_centre
	src_centre = np.vstack((int(num[5:]),ra_centre,dec_centre,start_time,stop_time))
	np.savetxt("%s/%s/ME/centre.dat"%(Wpath,num),src_centre.T,fmt="%s")
	#####################run the script of hxmtsoft###################################################################
	print num
	if not os.path.exists("%s/%s/ME/me_pi.fits"%(Wpath,num)):
		os.system('mepical evtfile=%s tempfile=%s outfile=%s/%s/ME/me_pi.fits'%(meevtfile,metempfile,Wpath,num))
		print num," : pi end"
	if not os.path.exists("%s/%s/ME/me_dt.fits"%(Wpath,num)):
		os.system('megrade evtfile=%s/%s/ME/me_pi.fits outfile=%s/%s/ME/me_grade.fits deadfile=%s/%s/ME/me_dt.fits binsize=1'%(Wpath,num,Wpath,num,Wpath,num))
		print num," : grade end"
	if not os.path.exists("%s/%s/ME/me_gti.fits"%(Wpath,num)):
		if flag5!=1:
			os.system('hxmtehkgen orbfile=%s attfile=%s outfile=%s/%s/ehk.fits step_sec=0.25 leapfile=/home/hxmt/guanj/zhaohsV2/hxmtehkgen/refdata/leapsec.fits rigidity=/home/hxmt/guanj/zhaohsV2/hxmtehkgen/refdata/rigidity_20060421.fits saafile=/home/hxmt/guanj/zhaohsV2/hxmtehkgen/SAA/SAA.fits'%(orbitfile,attfile,Wpath,num))
			os.system('megtigen tempfile=%s ehkfile=%s/%s/ehk.fits outfile=%s/%s/ME/gtiv2.fits defaultexpr=NONE expr="ELV>5&&COR3>=8&&T_SAA>=200&&TN_SAA>=100&&SAA_FLAG==0&&SUN_ANG>=10&&MOON_ANG>=5&&ANG_DIST<=359&&(SAT_LAT<31||SAT_LAT>38)&&(SAT_LON>245||SAT_LON<228)&&(SAT_LAT>=-36.5&&SAT_LAT<=36.5)"'%(metempfile,Wpath,num,Wpath,num)) 	
		else:
			os.system("cp %s %s/%s/ehk.fits"%(ehkfile,Wpath,num))
			os.system('megtigen tempfile=%s ehkfile=%s/%s/ehk.fits outfile=%s/%s/ME/gtiv2.fits defaultexpr=NONE expr="ELV>5&&COR>=8&&T_SAA>=200&&TN_SAA>=100&&SAA_FLAG==0&&SUN_ANG>=10&&MOON_ANG>=5&&ANG_DIST<=359&&(SAT_LAT<31||SAT_LAT>38)&&(SAT_LON>245||SAT_LON<228)&&(SAT_LAT>=-36.5&&SAT_LAT<=36.5)"'%(metempfile,Wpath,num,Wpath,num)) 
		print num," : o-gti & ehk end"
		regti_v2.regti(meevtfile,'%s/%s/ME/'%(Wpath,num))
		print num," : r-gti end"
	if not os.path.exists("%s/%s/bd_%s.fits"%(Wpath,num,num)):
		mescreen_newbd.mescreen_newbd('%s'%Wpath,'%s'%num)
		print num," : bd end"
	if not os.path.exists("%s/%s/ME/me_screen.fits"%(Wpath,num)):
		os.system('mescreen evtfile=%s/%s/ME/me_grade.fits gtifile=%s/%s/ME/me_gti.fits outfile=%s/%s/ME/me_screen.fits baddetfile=%s/%s/bd_%s.fits userdetid="0-53"'%(Wpath,num,Wpath,num,Wpath,num,Wpath,num,num))
		print num," : screen end"
	#os.system('mescreen evtfile=%s/%s/ME/me_grade.fits gtifile=%s/%s/ME/me_gti.fits outfile=%s/%s/ME/me_screen.fits baddetfile=${HEADAS}/refdata/medetectorstatus.fits userdetid="0-53"'%(Wpath,num,Wpath,num,Wpath,num))
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
	melcgen evtfile=%s/%s/ME/me_screen.fits deadfile=%s/%s/ME/me_dt.fits outfile=%s/%s/ME/me_${box}_small userdetid="${small}" starttime=0 stoptime=0 minPI=68 maxPI=631 binsize=1 deadcorr=yes
	melcgen evtfile=%s/%s/ME/me_screen.fits deadfile=%s/%s/ME/me_dt.fits outfile=%s/%s/ME/me_${box}_blind userdetid="${blind}" starttime=0 stoptime=0 minPI=68 maxPI=631 binsize=1 deadcorr=yes
	done'''%(Wpath,num,Wpath,num,Wpath,num,Wpath,num,Wpath,num,Wpath,num)
	#melcgen evtfile=%s/%s/ME/me_screen.fits deadfile=%s/%s/ME/me_dt.fits outfile=%s/%s/ME/me_${box}_big userdetid="${big}" starttime=0 stoptime=0 minPI=68 maxPI=631 binsize=1 deadcorr=yes
	if not os.path.exists("%s/%s/ME/me_2_blind_g0_46.lc"%(Wpath,num)):
		print me_lc_cmd
		#print num
		os.system(me_lc_cmd)
		print num," : o-lc end"
	os.system('rm -r %s/pfiles_%s'%(Ppath,num))
	print num," 1.Soft Cal: End"
	###########################
	if not os.path.exists("%s/%s/ME/me_2_blind.lc"%(Wpath,num)):
		errorCorrect.errorCorrect('%s'%Wpath,'%s'%num)
		print num," 2.Correct LC: Error Correction."
	if not os.path.exists("%s/%s/ME/me_lc_box2_small_bkg.fits"%(Wpath,num)):
		nbkghesmtpoly2.poly_bkg('%s'%Wpath,'%s'%Wpath,'%s'%num)
		print num, "3.BKG Cal: LC without BKG and BKG value."
	if not os.path.exists("%s/%s/ME/me_lc_box2_small.fits"%(Wpath,num)):
		gen_lc_me.gen_lc_me("%s/%s/ehk.fits"%(Wpath,num),"%s/%s/ME/me_gti.fits"%(Wpath,num),'%s/%s/ME/'%(Wpath,num),"%s"%num)
		print num, "4.Auto LC Cut 1: Manual Preparation."
	if not os.path.exists("%s/%s/ME/me_lc_box2_small_cut.fits"%(Wpath,num)):
		check_reatt_auto.reatt('%s/%s/att.fits'%(Wpath,num),'%s/%s/ME/'%(Wpath,num))
		print num, "5.Auto LC Cut 2: Auto Cut."


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

	###Wpath='/sharefs/hbkg/user/saina/data295/P0101295/'
	###path='/hxmt/work/HXMT-DATA/1L/A01/P0101295/'
	#lclist= os.listdir(Wpath)
	#lclist.sort()
	###lclist='P010129500105'
	merun_v2(path,Wpath,ObsID)

