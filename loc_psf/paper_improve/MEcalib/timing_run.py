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
sys.path.append("/sharefs/hbkg/user/luoqi/HXMT_SCAN/ME/ME_SCAN/lcfit")
from write_007xml_index import *

import subprocess###
def exec_cmd(cmd,my_env):
    """Run shell command"""
    p = subprocess.Popen(cmd, stdin=subprocess.PIPE,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT,
                         env=my_env,
                         shell=True)

    log_content = p.communicate()[0]
    #print(log_content.decode('utf-8'))

    return p.returncode, log_content

###os.system("source /hxmt/home/hxmtsoft2/hxmtsoft_v2.sh")
def mkdir_try(dirname):
	if os.path.exists(dirname) == 0:
		try:
			os.makedirs(dirname)
		except OSError:
			print('Wrong with make dir:\n'+dirname)

def merun_v2(path, Wpath, ObsID,Erange):
	if Erange=='7_12':
		program_tree = '/sharefs/hbkg/user/luoqi/psfl/calib/7_12/genlc'
		minpi = 68  ###energy =7 keV
		maxpi = 154  ###energy = 12 keV
	elif Erange=='7_20':
		program_tree = '/sharefs/hbkg/user/luoqi/psfl/calib/7_20/genlc'
		minpi = 68  ###energy =7 keV
		maxpi = 290  ###energy = 20 keV
	elif Erange=='7_40':
		program_tree = '/sharefs/hbkg/user/luoqi/psfl/calib/7_20/genlc'
		minpi = 68
		maxpi = 631
	else:
		print("no exist input Erange, please change code yourself")
		sys.exit(0)
	###--------------------mkdirs---------------------
	OrgWpath = Wpath + '/Org'  ###
	NetWpath = Wpath + '/Net'  ###
	mkdir_try(OrgWpath)
	mkdir_try(NetWpath)
	Org_obspath = "%s/%s"%(OrgWpath, ObsID)###
	Net_obspath = "%s/%s"%(NetWpath, ObsID)###
	Org_attpath = "%s/%s"%(OrgWpath, 'Att')###
	Org_gtipath = "%s/%s"%(OrgWpath, 'GTI')###
	Org_ehkpath = "%s/%s"%(OrgWpath, 'EHK')###
	Org_pipath = "%s/%s"%(OrgWpath, 'PI')###
	Org_screenpath = "%s/%s"%(OrgWpath, 'Screen')###
	mkdir_try(Org_obspath)
	mkdir_try(Net_obspath)
	mkdir_try(Org_attpath)
	mkdir_try(Org_gtipath)
	mkdir_try(Org_ehkpath)
	mkdir_try(Org_pipath)
	mkdir_try(Org_screenpath)
	###--------------------mkdirs---------------------
	###--------------------outfiles-------------------
	Wattfile = "%s/%s_Att.fits" % (Org_attpath, ObsID)  ###perfect
	Wgtifile = "%s/%s_gtiv2.fits" % (Org_gtipath, ObsID)  ###
	centre_txtfile = "%s/%s/centre.dat" % (OrgWpath, ObsID)  ###
	Wpifile = "%s/%s_me_pi.fits" % (Org_pipath, ObsID)  ###perfect
	Wdeadfile = "%s/%s_me_dead.fits" % (Org_obspath, ObsID)  ###dead time file
	Wgradefile = "%s/%s_me_grade.fits" % (Org_obspath, ObsID)  ###
	Wehkfile = "%s/%s_ehk.fits" % (Org_ehkpath, ObsID)  ###perfect
	Wme_gtifile = "%s/%s_me_gti.fits" % (Org_gtipath, ObsID)  ###perfect
	Wbaddetfile = "%s/%s_me_bdet.fits" % (Org_obspath, ObsID)  ###bad det
	Wscreenfile = "%s/%s_me_screen.fits" % (Org_screenpath, ObsID)  ###
	###--------------------outfiles-------------------

	stObsID = str(ObsID[:-2])###short ObsID
	#Wpath = "/hxmt/work/HXMT-DATA/1L/A02/P0211007/"+stObsID+"/"
	Npath = path[:22]+'N/'+path[28:]
	#Wpath = "/sharefs/hbkg/user/saina/data294/P0211007/"
	print path,Wpath
	###Change Ppath by LQ
	###Ppath='/sharefs/hbkg/user/saina/pfiles_small/'
	Ppath = '/sharefs/hbkg/user/luoqi/psfl/calib/pfiles_small'
	if os.path.exists(Ppath +"/pfiles_" + ObsID):
		os.system("rm -r %s/pfiles_%s" % (Ppath, ObsID))
	my_env = os.environ
	os.system("mkdir %s/pfiles_%s" % (Ppath, ObsID))
	my_env['PFILES']="%s/pfiles_%s;/home/hxmt/hxmtsoft2/soft/install/x86_64-unknown-linux-gnu-libc2.12/syspfiles"%(Ppath, ObsID)
	os.system("echo $PFILES")
	######################get all the files and folders##########################################################
	attfile_list = glob.iglob(r'%s/%s/*/*_Att_*'%(path,stObsID))
	meevtfile_list = glob.glob(r'%s/%s/%s*/*/*ME-Evt*' % (path, stObsID, ObsID))
	metempfile_list = glob.glob(r'%s/%s/%s*/*/*ME-TH*' % (path, stObsID, ObsID))
	orbitfile_list = glob.glob(r'%s/%s/*/*_Orbit_*'%(path,stObsID))
	ehk1N_path = "/sharefs/hbkg/user/cwang/ehk"
	ehkfile_list1 = glob.iglob(r'%s/ehk%s*'%(ehk1N_path,stObsID[2:]))##/sharefs/hbkg/user/cwang/ehk stObsID[2:]#glob.iglob(r'%s/%s/*/*EHK*'%(path,stObsID))
	ehkfile_list2 = glob.iglob(r'%s/%s/*/*EHK*'%(path,stObsID))
	print(ehkfile_list1)
	print(ehkfile_list2)
	nattfile_list = glob.iglob(r'%s/%s/*/*_Att_*'%(Npath,stObsID))
	norbitfile_list = glob.glob(r'%s/%s/*/*_Orbit_*'%(Npath,stObsID))
	print "%s: files loaded"%stObsID
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
	for i in ehkfile_list1:
		if re.findall(r"VU",i)==[]:
			ehkfile1 = i
			flag5=1
	for i in ehkfile_list2:
		if re.findall(r"VU",i)==[]:
			ehkfile2 = i
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
		print "%s %s: 1st -- Files error"%(stObsID, ObsID)

	###---old dir code position

	try:
		print orbitfile
	except NameError:
		print "%s %s can not be calculated: orbit"%(stObsID, ObsID)
		if (nflag4==1)&(nflag1==1):
			print "%s %s 1N orbit test"%(stObsID, ObsID)
			attfile=nattfile
			orbitfile = "%s/%s/%s_Orbit_1N.fits" % (OrgWpath, ObsID, ObsID)###
			fix_eci_1N.fix_eci(norbitfile,orbitfile)
		else:
			return 0
	try:
		print attfile
	except NameError:
		print "%s can not be calculated: att" % ObsID
		return 0
	try:
		print meevtfile
	except NameError:
		print "%s can not be calculated: evt" % ObsID
		return 0
	try:
		print metempfile
	except NameError:
		print "%s can not be calculated: th" % ObsID
		return 0
	try:
		print ehkfile1
		ehkfile = ehkfile1
	except NameError:
		try:
			print ehkfile2
			ehkfile = ehkfile2
		except NameError:
			print "No ehk file, need calculate by soft"
			flag5=0
	print "---------------------------------------------------------------------------------------"
	###################copy files to the floders##############################################################
	os.system("cp %s %s"%(attfile,Wattfile))
	os.system("rm %s"%(Wgtifile))
	evt_im = pf.open("%s"%(meevtfile))
	ename = evt_im
	ra_centre = ename[1].header['RA_PNT']
	dec_centre = ename[1].header['DEC_PNT']
	start_time = ename[1].header['TSTART']
	stop_time = ename[1].header['TSTOP']
	print "centre: ",ra_centre,dec_centre
	src_centre = np.vstack((int(ObsID[5:]), ra_centre, dec_centre, start_time, stop_time))
	np.savetxt(centre_txtfile, src_centre.T, fmt="%s")
	#####################run the script of hxmtsoft###################################################################
	print ObsID
	if 1:###not os.path.exists(Wpifile):
		cmd = 'mepical evtfile=%s tempfile=%s outfile=%s' % (meevtfile, metempfile, Wpifile)
		print(cmd)
		os.system(cmd)
		print ObsID, " : pi end"
	if 1:###not os.path.exists(Wdeadfile):
		cmd = 'megrade evtfile=%s outfile=%s deadfile=%s binsize=1' % (Wpifile, Wgradefile, Wdeadfile)
		os.system(cmd)
		print ObsID, " : grade end"

	os.system("cp %s %s" % (ehkfile, Wehkfile))

	###os.system('rm %s' % (Wme_gtifile))
	if 1:###not os.path.exists(Wgtifile):###----maybe need flag5, hxmtehkgen
		#if flag5!=1:
			#os.system("cp %s %s/%s/ehk.fits"%(ehkfile,Wpath,ObsID))#######after SAA_FLAG###SUN_ANG>=10&&MOON_ANG>=5&&
		cmd = 'megtigen tempfile=%s ehkfile=%s outfile=%s defaultexpr=NONE expr="ELV>5&&COR>=8&&T_SAA>=200&&TN_SAA>=100&&SAA_FLAG==0&&SUN_ANG>=10&&MOON_ANG>=5&&ANG_DIST<=359&&(SAT_LAT<31||SAT_LAT>38)&&(SAT_LON>245||SAT_LON<228)&&(SAT_LAT>=-36.5&&SAT_LAT<=36.5)"' % (
			metempfile, Wehkfile, Wgtifile)
		###cmd = 'megtigen tempfile=%s ehkfile=%s outfile=%s defaultexpr=YES'%(
		###	metempfile, Wehkfile, Wgtifile)
		print(cmd)
		with open(program_tree + '/gti_manual.txt', 'a+') as f:
			print>>f,cmd
		cmd_tem = ';export HEADASNOQUERY=;export HEADASPROMPT=/dev/null;' + cmd
		exec_code, exec_log = exec_cmd(cmd_tem, my_env)
		###os.system(cmd)
			#os.system('hxmtehkgen orbfile=%s attfile=%s outfile=%s/%s/ehk.fits step_sec=0.25 leapfile=/home/hxmt/guanj/zhaohsV2/hxmtehkgen/refdata/leapsec.fits rigidity=/home/hxmt/guanj/zhaohsV2/hxmtehkgen/refdata/rigidity_20060421.fits saafile=/home/hxmt/guanj/zhaohsV2/hxmtehkgen/SAA/SAA.fits'%(orbitfile,attfile,Wpath,ObsID))
			#os.system('megtigen tempfile=%s ehkfile=%s/%s/ehk.fits outfile=%s/%s/ME/gtiv2.fits defaultexpr=NONE expr="ELV>5&&COR3>=8&&T_SAA>=200&&TN_SAA>=100&&SAA_FLAG==0&&SUN_ANG>=10&&MOON_ANG>=5&&ANG_DIST<=359&&(SAT_LAT<31||SAT_LAT>38)&&(SAT_LON>245||SAT_LON<228)&&(SAT_LAT>=-36.5&&SAT_LAT<=36.5)"'%(metempfile,Wpath,ObsID,Wpath,ObsID))
		#else:
			#os.system("cp %s %s/%s/ehk.fits"%(ehkfile,Wpath,ObsID))
			#os.system('megtigen tempfile=%s ehkfile=%s/%s/ehk.fits outfile=%s/%s/ME/gtiv2.fits defaultexpr=NONE expr="ELV>5&&COR>=8&&T_SAA>=200&&TN_SAA>=100&&SAA_FLAG==0&&SUN_ANG>=10&&MOON_ANG>=5&&ANG_DIST<=359&&(SAT_LAT<31||SAT_LAT>38)&&(SAT_LON>245||SAT_LON<228)&&(SAT_LAT>=-36.5&&SAT_LAT<=36.5)"'%(metempfile,Wpath,ObsID,Wpath,ObsID))
		print ObsID, " : o-gti & ehk end"
		regti_v2.regti(meevtfile, Wpath, ObsID)###calc total time
		print ObsID, " : r-gti end"

	if not os.path.exists(Wbaddetfile):
		mescreen_newbd.mescreen_newbd('%s' % Wpath,'%s' % ObsID)###
		print ObsID, " : bd end"

	if 1:###not os.path.exists(Wscreenfile):
		cmd = 'mescreen evtfile=%s gtifile=%s outfile=%s baddetfile=%s userdetid="0-53"' % (Wgradefile, Wme_gtifile, Wscreenfile, Wbaddetfile)
		print(cmd)
		os.system(cmd)
		print ObsID, " : screen end"
	#os.system('mescreen evtfile=%s/%s/ME/me_grade.fits gtifile=%s/%s/ME/me_gti.fits outfile=%s/%s/ME/me_screen.fits baddetfile=${HEADAS}/refdata/medetectorstatus.fits userdetid="0-53"'%(Wpath,ObsID,Wpath,ObsID,Wpath,ObsID))
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
	melcgen evtfile=%s deadfile=%s outfile=%s/%s/me_${box}_small userdetid="${small}" starttime=0 stoptime=0 minPI=%s maxPI=%s binsize=1 deadcorr=yes
	melcgen evtfile=%s deadfile=%s outfile=%s/%s/me_${box}_blind userdetid="${blind}" starttime=0 stoptime=0 minPI=%s maxPI=%s binsize=1 deadcorr=yes
	done'''%(Wscreenfile, Wdeadfile, NetWpath, ObsID, minpi, maxpi, Wscreenfile, Wdeadfile, NetWpath, ObsID, minpi, maxpi)###change outfile name----
	#melcgen evtfile=%s/%s/ME/me_screen.fits deadfile=%s/%s/ME/me_dt.fits outfile=%s/%s/ME/me_${box}_big userdetid="${big}" starttime=0 stoptime=0 minPI=68 maxPI=631 binsize=1 deadcorr=yes
	if 1:###not os.path.exists("%s/%s/me_2_blind_g0_46.lc"%(NetWpath, ObsID)):
		print me_lc_cmd
		#print ObsID
		os.system(me_lc_cmd)
		print ObsID, " : o-lc end"
	os.system('rm -r %s/pfiles_%s' % (Ppath, ObsID))
	print ObsID, " 1.Soft Cal: End"
	###########################
	if 1:###not os.path.exists("%s/%s/me_2_blind.lc"%(NetWpath, ObsID)):
		errorCorrect.errorCorrect('%s' % Wpath,'%s' % ObsID)
		print ObsID, " 2.Correct LC: Error Correction."
	if 1:###not os.path.exists("%s/%s/me_lc_box2_small_bkg.fits"%(NetWpath, ObsID)):
		nbkghesmtpoly2.poly_bkg('%s' % NetWpath,'%s' % Wpath,'%s' % ObsID)
		print ObsID, "2.BKG Cal: LC without BKG and BKG value."
	if 1:###not os.path.exists("%s/%s/me_lc_box2_small.fits"%(NetWpath, ObsID)):
		gen_lc_me.gen_lc_me(Wehkfile, Wme_gtifile, '%s' % (Wpath), "%s" % ObsID)
		print ObsID, "3.Auto LC Cut 1: Manual Preparation."
	if 1:###not os.path.exists("%s/%s/me_lc_box2_small_cut.fits"%(NetWpath, ObsID)):
		check_reatt_auto.reatt(Wattfile, Wpath, ObsID)
		print ObsID, "4.Auto LC Cut 2: Auto Cut."

	xml_file = write_xml(Wpath, ObsID)
	lcfit_command = 'source /sharefs/hbkg/user/luoqi/home/scan_MEfit_env.sh;sh %s/lcfit/sub_gps_lcfit.sh %s' % (program_tree,
		xml_file)
	###lcfit_command = 'source /sharefs/hbkg/user/luoqi/home/scan_MEfit_env.sh;python /sharefs/hbkg/user/luoqi/HXMT_SCAN/ME/ME_SCAN/lcfit/test_all_typeSub.py %s'%(xml_file)
	with open(program_tree + '/genlc%s_success.txt'%Erange, 'a+') as f:
		print>>f,ObsID
	with open(program_tree + '/lcfit%s_misson.sh'%Erange, 'a+') as f:
		print>>f,lcfit_command
	###os.system(lcfit_command)



if __name__ == "__main__":

	if len(sys.argv) < 2:
		ObsID = input("input keywords like P010129400101 :")
	elif len(sys.argv) >= 2:
		ObsID = sys.argv[1]

	Erange = sys.argv[2]
	scan_tree = '/sharefs/hbkg/data/SCAN/luoqi/calib/%s' % Erange
	Wpath = scan_tree  ###'/sharefs/hbkg/data/SCAN/ME'
	path = '/hxmt/work/HXMT-DATA/1L/A%s/%s/' % (ObsID[1:3], ObsID[:-5])


	###path = '/hxmt/work/HXMT-DATA/1L/A01/P0101295'
	###ObsID = 'P010129500103'

	###Wpath='/sharefs/hbkg/user/saina/data294/P0211007/'
	###Wpath='/hxmt/work/HXMT-DATA/1L/A02/P0211007/'
	#lclist= os.listdir(Wpath)
	#lclist.sort()

	###lclist='P021100715501'
	merun_v2(path,Wpath,ObsID,Erange)

