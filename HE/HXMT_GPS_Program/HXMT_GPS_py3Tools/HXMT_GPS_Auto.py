import os,sys
#sys.path.append("/sharefs/hbkg/user/nangyi/HXMT_GPS_Tools")
#sys.path.append("/Users/songxy/work/ihep4/HXMT_GPS_Tools")
sys.path.append("/sharefs/hbkg/user/luoqi/GRB/work/ihep4/HXMT_GPS_Tools")
import multiprocessing as mp
from Lctools_201903 import FindLC_basefile
import time

data1L_path = "/hxmt/work/HXMT-DATA/1L/A01/P0101295/"#"/hxmt/work/HXMT-DATA/1L/A02/P0201013/"
genlc_logfile = 'Genlc_process.log'
normlc_logtmp = 'Genlc_tmp.log'
normlc_logfile = 'Genlc_succes.log'
srpath = '/sharefs/hbkg/users/luoqi/GRB/work/ihep4/data/HE_GPS_SML/'#'/sharefs/hbkg/users/luoqi/GRB/work/ihep4/data/HE_GPS_SML/'

def fucgenlc(path):
    print str(path)
    temp = os.system("python GenLc.py %s"%(str(path)))
    return [path,temp]

def fuc_normfit(path):
    print str(path)
    pfile = "/home/hxmt/nangyi/pfiles/%s"%path
    pfilepath = "%s;/hxmt/soft/Astro/heasoft/x86_64-pc-linux-gnu-libc2.12/syspfiles"%pfile
    commond = 'rm /hxmt/work/HXMT_scan_data/src_list/HE/*%s*.txt;mkdir %s;cd %s%s;rm `ls %s%s/* | grep -v config_he.xml`;cp /sharefs/hbkg/user/nangyi/run/simpleknow.py %s%s/;source /hxmt/work/USERS/nangyi/envV2.sh ;export PFILES="%s"; python simpleknow.py config_he.xml;rm -rf %s;\n'%(path,pfile,srpath,path,srpath,path,srpath,path,pfilepath,pfile)
    print commond
    temp = os.system(commond)
    return [path,temp]

def fuc_grbfit(path):
    print str(path)
    pfile = "/home/hxmt/nangyi/pfiles/%s"%path
    pfilepath = "%s;/hxmt/soft/Astro/heasoft/x86_64-pc-linux-gnu-libc2.12/syspfiles"%pfile
    commond = 'rm /hxmt/work/HXMT_scan_data/src_list/HE/*%s*.txt;mkdir %s;cd %s%s;rm `ls %s%s/* | grep -v config_he.xml`;cp /sharefs/hbkg/user/nangyi/run/simplegrb.py %s%s/;source /hxmt/work/USERS/nangyi/envV2.sh ;export PFILES="%s"; python simplegrb.py config_he.xml;rm -rf %s;\n'%(path,pfile,srpath,path,srpath,path,srpath,path,pfilepath,pfile)
    print commond
    temp = os.system(commond)
    return [path,temp]

def Run_pool(fuc, lists):
    pool = mp.Pool(processes = 10)
    res = pool.map(fuc, lists)
    pool.close()
    pool.join()
    print "Pool Finished!"

def deal_newdata(lists):
    print 'Generate LC:::',lists
    Run_pool(fucgenlc,lists)
    with open(genlc_logfile,'a')as f:
        for i in lists:
            tmptm = str(time.strftime('%Y.%m.%d-%H:%M',time.localtime(time.time())))
            print>>f,i,'Processed!',tmptm
    norm_lc = []
    try:
        with open(normlc_logtmp,'r')as f:
            for i in f:
                norm_lc.append(i[:13])
        os.system('rm %s'%normlc_logtmp)
        with open(normlc_logfile,'a')as f:
            for i in norm_lc:
                emptm = str(time.strftime('%Y.%m.%d-%H:%M',time.localtime(time.time())))
                print>>f,i,' Normal-Generated!',tmptm
    except:
        pass
    print norm_lc
    norm_lclist = [i[5:13] for i in norm_lc]
    if len(norm_lc)>0:
        print 'Begin to Fit>>>'
        Run_pool(fuc_normfit,norm_lclist)


def run():
    #folder_now = os.listdir(data1L_path)
    genlc_loglist = []
    with open(genlc_logfile,'r')as f:
        for i in f:genlc_loglist.append(i[:11])
        print(i[:11],i)
    print('ok')
    '''
    mins = set(folder_now)-set(genlc_loglist)
    folder_acess = []
    for i in mins:
        print i
        try:
            temp = FindLC_basefile(i+'01',data1L_path+i+'/',INST='HE',EHKFILE=False)
            folder_acess.append(i+'01')
            print 'New folder accessible',i
        except:
            print 'New folder without enough FILES'
            pass
        #folder_acess.sort(reverse=True)
        #print folder_acess
    if len(folder_acess)>0:
        folder_acess.sort(reverse=True)
        print 'New Obervation Accessible!',folder_acess
        deal_newdata(folder_acess)
    else:
        print 'There is no New Obervation!'
        time.sleep(1800)
    '''



if __name__ == "__main__":
    while True:
        run()

