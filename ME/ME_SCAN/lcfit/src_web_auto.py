#!/hxmt/soft/Develop/anaconda2/bin/python
import numpy as np
import os

#list1=[]
#Wpath = '/sharefs/hbkg/user/saina/data294/'
#Wpath1='/sharefs/hbkg/user/saina/data294/P0101294/'
#Wpath2='/sharefs/hbkg/user/saina/data294/P0211007/'
#Wpath3='/sharefs/hbkg/user/saina/data294/P0201013/'
#dirlist = os.listdir(Wpath1)+os.listdir(Wpath2)+os.listdir(Wpath3)
#dirlist.sort()
def mkdir_try(dirname):
    if os.path.exists(dirname) == 0:
        try:
            os.makedirs(dirname)
        except OSError:
            print('Wrong with make dir:\n'+dirname)

def dataout(Wpath, ObsID):
	Program_dir = '/sharefs/hbkg/user/luoqi/HXMT_SCAN/ME/ME_SCAN'
	#path = Wpath
	Midd_obspath = Wpath + '/Midd/' + ObsID  ###
	fProd_obspath = Wpath + '/Prod/' + ObsID  ###
	#print Wpath,ObsID
	if not os.path.exists("%s/%s_data_3box.dat"%(Midd_obspath, ObsID)):
		print "%s: step1 not end" % ObsID
		return 0
	if not os.path.exists("%s/%s_source_ub.txt"%(Midd_obspath, ObsID)):###???
		print "%s: step2 not end" % ObsID
		return 0
	print "%s: start" % ObsID
	src_list2 = np.loadtxt("%s/%s_source_ub.txt" % (Midd_obspath, ObsID), delimiter='\t', dtype=bytes).astype(str)
	if os.path.exists("%s/%s_source_bright.txt"%(Midd_obspath, ObsID)):
		src_listo = np.loadtxt("%s/%s_source_bright.txt" % (Midd_obspath, ObsID), delimiter='\t', dtype=bytes).astype(str)
		src_listm = np.loadtxt("%s/%s_source_moreinfo.txt" % (Midd_obspath, ObsID), delimiter='\t', dtype=bytes).astype(str)
		if len(src_listm.shape)<=1:
			if src_listm[6].astype(float)>0:
				src_list1 = src_listo
			else:
				src_list1 = [0]
		else:
			src_list1 = src_listo[src_listm[:,7].astype(float)>0]
		#print src_list1
		if len(src_list1)>1:
			if len(src_list1.shape)<=1:
				src_list1=[src_list1]
			src_list=np.r_[src_list1,src_list2]
		else:
			src_list=src_list2
	else:
		src_list=src_list2
	i2=0
	src_list_new=src_list[np.lexsort(src_list[:,1::-2].astype(int).T)]
	for i2 in range(len(src_list_new)):
		#src_list_new[i2][0]=i2+1
		if (src_list_new[i2][7].astype(float)==99999999)|(src_list_new[i2][7]==None)|((src_list_new[i2][7].astype(float)==99999999)&(src_list_new[i2][5].astype(float)<=0))|((src_list_new[i2][7].astype(float)<=0)&(src_list_new[i2][5].astype(float)<=0)):
			#src_list_new[i][7]=99999999
			src_list_new=np.delete(src_list_new,i2,axis=0)
			i2=i2-1
		i2=i2+1
	if len(src_list_new)==0:
		print '%s: list error' % ObsID
		return 0
	num0=np.array(range(len(src_list_new[:,0])))+1
	src_list_new[:,0]=num0
	###f=open("%s/src_ME/src_%s.txt" % (Program_dir,ObsID), "w+")###Midd_obspath
	f = open("%s/src_%s.txt" % (Midd_obspath, ObsID), "w+")
	for num in src_list_new:
		filestr1='\t'.join(num)
		f.writelines(filestr1+'\n')
	f.close()
	mkdir_try(fProd_obspath)
	os.system('cp %s/src_%s.txt %s/src_%s.txt'% (Midd_obspath, ObsID, fProd_obspath, ObsID))
	###after should copy#os.system('cp %s/src_ME/src_%s.txt /hxmt/work/HXMT_scan_data/src_list/ME/' % (Program_dir,ObsID))

if __name__ == "__main__":
	lclist='P020101300201'
	Wpath='/sharefs/hbkg/user/saina/data294/'+lclist[:-5]+'/'+lclist+'/ME/result/'
	#lclist= os.listdir(Wpath)
	#lclist.sort()
	dataout(Wpath,lclist)

