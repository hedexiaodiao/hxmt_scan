import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits as pf
import sys,os
from xml.etree.ElementTree import ElementTree,Element
import time
Program_dir = '/sharefs/hbkg/user/luoqi/HXMT_SCAN/ME/ME_SCAN'
def mkdir_try(dirname):
    if os.path.exists(dirname) == 0:
        try:
            os.makedirs(dirname)
        except OSError:
            print('Wrong with make dir:\n'+dirname)

def write_xml(Wpath, ObsID):
	OrgWpath = Wpath + '/Org'  ###
	NetWpath = Wpath + '/Net'  ###
	Org_obspath = "%s/%s" % (OrgWpath, ObsID)  ###
	Net_obspath = "%s/%s" % (NetWpath, ObsID)  ###
	Midd_obspath = Wpath + '/Midd/' + ObsID###
	mkdir_try(Midd_obspath)
	Org_attpath = "%s/%s" % (OrgWpath, 'Att')  ###
	#Wattfile = "%s/%s_Att.fits" % (Org_attpath, ObsID)###
	Wcut_attfile = "%s/%s_cut_Att.fits" % (Org_attpath, ObsID)
	os.system('cp %s/centre.dat %s/'%(Org_obspath,Midd_obspath))
	if 1:###os.path.exists("%s/"%Net_obspath):
		xmlpath = "%s/lcfit/config_me_small.xml"%Program_dir###
		filestr = '%s/'%(Midd_obspath)
		#path1="/hxmt/work/HXMT-DATA/1L/A02/P0211007/P0211007"+ObsID+"/"
		#path2="/sharefs/hbkg/data/DATA/1L/A01/Temp_P0101294/%s/O010129%s"%(day,ObsID)
		bkgpath="%s"%Net_obspath
		tree = ElementTree()
		tree.parse(xmlpath) 
		###innodes = tree.findall("PATH/inpath")
		pathnodes = tree.findall("PATH/outpath")
		outnodes 		= tree.findall("Outfile_Name/outfilename")
		filenodes = tree.findall("Infile_Name/infilename/infilenamelist")
		#fileinpath = tree.findall("PATH/inpath/inpathnamelist")
		outfile = filestr + "config_me_small.xml"
		print outfile
		filenodes[0].text ='\n  '+ Net_obspath +'/me_lc_box0_small_cut.fits \n'###+ Net_obspath +###../../Net/'+ObsID+
		filenodes[1].text ='\n  '+ Net_obspath +'/me_lc_box1_small_cut.fits \n'###
		filenodes[2].text ='\n  '+ Net_obspath +'/me_lc_box2_small_cut.fits \n'###
		#filenodes[3].text ='\n  result/me_lc_box0_small_res.fits \n'
		#filenodes[4].text ='\n  result/me_lc_box1_small_res.fits \n'
		#filenodes[5].text ='\n  result/me_lc_box2_small_res.fits \n'
		filenodes[3].text ='\n  '+Wcut_attfile+' \n'
		#fileinpath[0].text = '\n ' + path1 +' \n'
		#fileinpath[1].text = '\n ' + path2 +' \n'
		#fileinpath[1].text = '\n ' + ObsID + ' \n'
		pathnodes[0].text ='\n  '+filestr+'  \n'###filestr
		###innodes[0].text = '\n  '+filestr+'  \n'
		tree.write(outfile, encoding="utf-8",xml_declaration=True)
	return outfile

if __name__ == "__main__":
	lclist = os.listdir("/sharefs/hbkg/data/SCAN/luoqi/ME/Net")
	###lclist=os.listdir("/sharefs/hbkg/user/saina/data294/P0201013/")
	lclist.sort()
	#a=[401025,401026,402023,403028,404029,405027,406020,406021,412022,413020,414027,415031,415032,416020,416021]
	#b=[12022,13020,14027,15031,15032,16020,16021]
	#lclist=["P0101294841"]
	Wpath = '/sharefs/hbkg/data/SCAN/luoqi/ME'
	###path= '/sharefs/hbkg/user/saina/data294/P0201013/'
	for i in lclist:
		write_xml(Wpath, i)
		print '%s is ok'%i
