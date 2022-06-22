import datetime
import time
import os
import glob
import numpy as np
from multiprocessing import Process
import sys
import smtplib
from email.mime.text import MIMEText
sys.path.append("/sharefs/hbkg/user/luoqi/HXMT_SCAN/ME/ME_SCAN/genlc")
###sys.Wpath.append("/sharefs/hbkg/user/saina/tar/P0201013")
import timing_run
import write_007xml_index
import threading

sem = threading.Semaphore(7) ###limit the maximum number of threads to 4

def fuc(path,fpath,lcpath):
	print str(path)
	#if not os.Wpath.exists("/sharefs/hbkg/user/saina/data294/P0211007/%s"%(str(Wpath))):
	if os.path.exists("%s/%s"%(fpath,str(path[:-2]))):
		print fpath,str(path), ": New calculation"
		Wpath=lcpath[:-4]
		write_007xml_index.write_xml(Wpath,str(path))
		timing_run.merun_v2(fpath,Wpath,str(path))
	else:
		print fpath,str(path),": No data"
	#else:
		#print str(Wpath),": Already done"


def doSth(list5,fpath,lcpath):
	plog = []
	for k in list5:
		for i in k:
			p = Process(target = fuc, args = (i,fpath,lcpath))
			p.start()
			plog.append(p)
	
		for p in plog:
			p.join()
			
	'''
	for k in list5:
		for i in k:
			fuc(i)
	'''
	print "Check END"
	time.sleep(60)



def send_mail(username, passwd, recv, title, content, mail_host='mail.ihep.ac.cn', port=25):
	msg = MIMEText(content,_subtype='plain',_charset='gb2312')###
	me = "<"+username+"@139.com>"
	msg['Subject'] = title
	msg['From'] = me
	msg['To'] = recv
	try:
		server = smtplib.SMTP()
		server.connect(mail_host)
		server.login(username, passwd)
		server.sendmail(me,recv,msg.as_string())
		server.close()
		return True
		print('email send success.')
	except :
		return False
	###smtp = smtplib.SMTP(mail_host, port=port)
	###smtp.login(username, passwd)
	###smtp.sendmail(username, recv, msg.as_string())
	###smtp.quit()



email_user = 'luoqi_ba'
email_pwd = 'luoqi@172'
maillist = 'liaojinyuan@ihep.ac.cn'
###maillist = 'luoqi@ihep.ac.cn'
mail_host = 'smtp.139.com'


def main(h=0, m=0):
	while True:
		stime=time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time()))
		#obin='P0101294'
		if len(sys.argv)<2:
			print "Need the ffn!"
			obin='P0'+str(input("ff ObsID (211007 211009 201013 101294 101295): "))
			print obin
			#sys.exit(1)
		else:
			print sys.argv[1]
			obin = sys.argv[1]
		fpath = "/hxmt/work/HXMT-DATA/1L/A%s/%s" % (obin[1:3], obin)
		'''
		obpath='/hxmt/work/HXMT-DATA/1L/'
		if obin[-1]=='7':
			fpath=obpath+'A02/'+obin###+'/'
		if obin[-1]=='4':
			fpath=obpath+'A01/'+obin###+'/'
		if obin[-1]=='9':
			fpath=obpath+'A02/'+obin###+'/'
		if obin[-1]=='5':
			fpath=obpath+'A01/'+obin###+'/'
		if obin[-1]=='3':
			fpath=obpath+'A02/'+obin###+'/'
		'''
		ldata_olist = glob.glob(r'%s/%s*/%s*/'%(fpath,obin,obin))
		ldatalist=[]
		for i in ldata_olist:
			ldatalist=np.append(ldatalist,i[-29:-16])
		
		###----------------------------------position for genlc------------------------------------
		'''
		if obin[-1]=='5':
			lcpath='/sharefs/hbkg/user/saina/data295/'+obin+'/'
		else:
			lcpath='/sharefs/hbkg/user/saina/data294/'+obin+'/'
		'''
		scan_path = '/sharefs/hbkg/data/SCAN/ME'
		lcpath = '/sharefs/hbkg/data/SCAN/ME/Net'

		lclist=os.listdir(lcpath)
		#ldatalist=os.listdir("/hxmt/work/HXMT-DATA/1L/A02/P0211007/")
		lclistre=[]
		for i in lclist:
			if not os.path.exists(lcpath+"/%s/me_lc_box2_small_cut.fits"%(str(i))):
				lclistre=np.append(lclistre,i)
		prodlist = []
		for j in ldatalist:
			prodfits_flag = os.path.exists(scan_path+"/Prod/%s/%s_FITDATA_ME.fits"%(str(j),str(j)))
			prodtxt_flag = os.path.exists(scan_path+"/Prod/%s/src_%s.txt"%(str(j),str(j)))
			if  (prodfits_flag and prodtxt_flag):
				prodlist = np.append(prodlist,j)

		
		#ldatalist=ldatalist[:700]
		# print("ldatalist",ldatalist)
		# mark=np.in1d(ldatalist,lclist,invert=True)
		# print("lclist",lclist)
		# print("mark",mark)
		# ldatalist=np.array(ldatalist)
		mins = list(set(ldatalist)-set(prodlist))
		ndatalist=mins###np.r_[mins,lclistre]
		ndatalist.sort()
		print("ndatalist",ndatalist)
		list5=[[]]*int(1+len(ndatalist)/5)
		for i in range(len(list5)):
			list5[i]=ndatalist[i*5:(i+1)*5]
			#print list5[i]
		
		rlist=[]
		doSth(list5,fpath,lcpath)
		for k in list5:
			for num in k:
				if os.path.exists(lcpath+"/%s/me_lc_box2_small_cut.fits"%(str(num))):
					rlist=np.append(rlist,num)
		
		#rlist = os.listdir("/sharefs/hbkg/user/saina/data294/P0211007/")
		#mark2=np.in1d(rlist,lclist,invert=True)
		#aplist=[]
		#if True in mark2:
			#aplist=np.array(rlist)[mark2]
		aplist=rlist
		#stime=time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time()))
		title = 'New data list: %s'%stime
		if len(aplist)>0:
			content = '%s'%aplist
			send_mail(email_user, email_pwd, maillist, title, content,mail_host)
		else:
			content = 'No new data today'
		###change apfile by LQ
		###apfile = '/sharefs/hbkg/user/saina/tar/'+obin+'/aplist.txt'
		try:
			os.makedirs('/sharefs/hbkg/user/luoqi/HXMT_SCAN/ME/apfiles/'+obin)
		except :
			pass
		apfile = '/sharefs/hbkg/user/luoqi/HXMT_SCAN/ME/apfiles/'+obin+'/aplist.txt'
		with open(apfile,'a') as f:
			for i in aplist:
				f.write('%s'%i+'\n')
		while True:
			now = datetime.datetime.now()
			#if (now.hour%3)==0 and (now.minute%5)==m:
			if (now.minute%5)==m:
				break
			time.sleep(1800)
		
with sem:
	main(0,0)
