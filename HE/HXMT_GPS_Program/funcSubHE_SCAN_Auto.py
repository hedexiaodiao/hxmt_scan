import os,sys
#sys.ObsID.append("/sharefs/hbkg/user/nangyi/HXMT_GPS_Tools")
#sys.ObsID.append("/Users/songxy/work/ihep4/HXMT_GPS_Tools")
###sys.ObsID.append("/sharefs/hbkg/user/luoqi/GRB/work/ihep4/HXMT_GPS_Tools")
sys.path.append(os.path.abspath(os.path.dirname(__file__))+'/HXMT_GPS_py3Tools')
import multiprocessing as mp
from Lctools_201903 import FindLC_basefile
import time
from exec_cmd import *



def fucgenlc(ObsID):
    print('==================>start genlc of obs:',str(ObsID))
    my_env = os.environ
    cmd_tem = 'source /sharefs/hbkg/user/luoqi/home/mypython;export HEADASNOQUERY=;export HEADASPROMPT=/dev/null;' +\
              "sh %s/HXMT_GPS_Program/sub_gps_allpro.sh %s"%(mkdir_path, ObsID)
    exec_code, exec_log = exec_cmd(cmd_tem, my_env)
    print(cmd_tem)
    print('##################################genlc ObsID, exec_code:', str(ObsID), exec_code)
    if exec_code==0:
        with open(genlc_logfile,'a')as f:
            tmptm = str(time.strftime('%Y.%m.%d-%H:%M',time.localtime(time.time())))
            print(str(ObsID),'Processed!',tmptm,file=f)
    return [ObsID,exec_code]


def fuc_grbfit(ObsID):
    print(str(ObsID))
    srpath = scan_dir + '/HE/GRB/Midd/'
    pfile = "%s/pfiles/%s"%(mkdir_path, ObsID)
    pfilepath = "%s;/hxmt/soft/Astro/heasoft/x86_64-pc-linux-gnu-libc2.12/syspfiles"%pfile
    commond = 'rm %s/src_list/*%s*.txt;mkdir %s;cd %s%s;cp %s/HXMT_GPS_Program/Lcfit.py %s%s/;source /sharefs/hbkg/user/luoqi/home/mypython;export PFILES="%s"; python %s/HXMT_GPS_Program/Lcfit.py %s%s/config_he.xml;rm -rf %s;\n' % (
        mkdir_path, ObsID, pfile, srpath, ObsID, mkdir_path, srpath, ObsID, pfilepath, mkdir_path, srpath, ObsID,
        pfile)
    ###commond = 'rm /hxmt/work/HXMT_scan_data/src_list/HE/*%s*.txt;mkdir %s;cd %s%s;rm `ls %s%s/* | grep -v config_he.xml`;cp /sharefs/hbkg/user/nangyi/run/simplegrb.py %s%s/;source /hxmt/work/USERS/nangyi/envV2.sh ;export PFILES="%s"; python simplegrb.py config_he.xml;rm -rf %s;\n'%(ObsID,pfile,srpath,ObsID,srpath,ObsID,srpath,ObsID,pfilepath,pfile)
    print(commond)
    temp = os.system(commond)
    return [ObsID, temp]

def Run_pool(fuc, lists):
    pool = mp.Pool(processes = 10)
    print()
    res = pool.map(fuc, lists)
    pool.close()
    pool.join()
    print("Pool Finished!")

def deal_newdata(lists):
    print('Generate LC:::',lists)
    Run_pool(fucgenlc,lists)
    with open(genlc_runningfile,'a')as f:
        for i in lists:
            tmptm = str(time.strftime('%Y.%m.%d-%H:%M',time.localtime(time.time())))
            print(i,'Processed!',tmptm,file=f)
    norm_lc = []
    try:
        with open(normlc_logtmp,'a')as f:###GRB mode will not execute
            for i in f:
                norm_lc.append(i[:13])
        os.system('rm %s'%normlc_logtmp)
    except:
        pass
    try:
        with open(normlc_logfile,'a')as f:
            for i in norm_lc:
                emptm = str(time.strftime('%Y.%m.%d-%H:%M',time.localtime(time.time())))
                print(i,' Normal-Generated!',tmptm,file=f)
    except:
        pass
    print(norm_lc)
    ##Run_pool(fuc_grbfit, lists)###GRB mode data for test

def run():
    circ_time1 = 0###2400
    circ_time2 = 1800
    exist_datalist = os.listdir(exist_proddir)
    exist_prodlist = [i[:13] for i in exist_datalist]
    folder_now = os.listdir(data1L_path)
    genlc_loglist = []
    with open(genlc_logfile,'r')as f:
        for i in f:genlc_loglist.append(i[:11])
    genlc_dealinglist = []###11
    with open(genlc_runningfile,'r')as f:
        for i in f:genlc_dealinglist.append(i[:11])
    mins1 = list(set(folder_now)-set(genlc_loglist))
    mins = list(set(mins1)-set(genlc_dealinglist))
    ##mins = mins[245:246]###
    print(mins)

    print('Number of ObsIDs:',len(mins))
    if len(mins)==0:
        print('There is no New Obervation!')
        time.sleep(circ_time2)###

    tog_num = 5
    steps = int(len(mins)/tog_num)
    last_IDs = len(mins)%tog_num
    for k in range(steps):
        folder_acess = []
        for i in mins[(tog_num*k):(tog_num*(k+1))]:
            print(i)
            for j in range(8):
                try:
                    temp = FindLC_basefile(i + '{:02d}'.format(j), data1L_path + i + '/', INST='HE', EHKFILE=False)
                    folder_acess.append(i + '{:02d}'.format(j))
                    print('New folder accessible', i + '{:02d}'.format(j))
                except:
                    print('New folder without enough FILES')
                    pass
            #folder_acess.sort(reverse=True)
            #print folder_acess
        if len(folder_acess)>0:
            folder_acess.sort(reverse=True)
            folder_acess0 = list(set(folder_acess)-set(exist_prodlist))
            print('New Obervation Accessible!',folder_acess0)
            deal_newdata(folder_acess0)
            time.sleep(circ_time1)
        else:
            print('There is no New Obervation!')
            time.sleep(circ_time2)###
    if last_IDs>0:
        folder_acess = []
        for i in mins[(tog_num*steps):len(mins)]:
            print(i)
            for j in range(8):
                try:
                    temp = FindLC_basefile(i + '{:02d}'.format(j), data1L_path + i + '/', INST='HE', EHKFILE=False)
                    folder_acess.append(i + '{:02d}'.format(j))
                    print('New folder accessible', i + '{:02d}'.format(j))
                except:
                    print('New folder without enough FILES')
                    pass
            #folder_acess.sort(reverse=True)
            #print folder_acess
        if len(folder_acess)>0:
            folder_acess.sort(reverse=True)
            folder_acess0 = list(set(folder_acess)-set(exist_prodlist))
            print('New Obervation Accessible!',folder_acess0)
            deal_newdata(folder_acess0)
            time.sleep(circ_time1)
        else:
            print('There is no New Obervation!')
            time.sleep(circ_time2)###
    if len(mins)<5:
        time.sleep(circ_time2)



if __name__ == "__main__":
    dir_strlist = []
    with open('dir_config.txt', 'r', encoding='utf-8') as f:
        for line in f.readlines():
            line = line.strip()
            dir_strlist.append(line)
    print(dir_strlist)
    mkdir_path = dir_strlist[0]
    scan_dir = dir_strlist[1]
    exist_proddir = '/sharefs/hbkg/data/SCAN/GPS_Production/HE/25_100keV'
    ###data1L_path = "/hxmt/work/HXMT-DATA/1L/A01/P0101295/"  # "/hxmt/work/HXMT-DATA/1L/A02/P0201013/"
    ###data1L_path = "/hxmt/work/HXMT-DATA/1L/A02/P0201013/"
    data1L_path = "/hxmt/work/HXMT-DATA/1L/A03/P0301240/"
    genlc_runningfile = mkdir_path+'/HXMT_GPS_Program/'+'subGenlc_running.log'
    genlc_logfile = mkdir_path+'/HXMT_GPS_Program/'+'subGenlc_process.log'
    normlc_logtmp = mkdir_path+'/HXMT_GPS_Program/'+'Genlc_tmp.log'
    normlc_logfile = mkdir_path+'/HXMT_GPS_Program/'+'Genlc_succes.log'
    ###srpath = '/sharefs/hbkg/users/luoqi/GRB/work/ihep4/data/HE_GPS_SML/'
    srpath = scan_dir + '/HE/Midd/'
    while True:
        run()


