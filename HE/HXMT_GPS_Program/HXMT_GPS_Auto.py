import os,sys
#sys.path.append("/sharefs/hbkg/user/nangyi/HXMT_GPS_Tools")
#sys.path.append("/Users/songxy/work/ihep4/HXMT_GPS_Tools")
###sys.path.append("/sharefs/hbkg/user/luoqi/GRB/work/ihep4/HXMT_GPS_Tools")
sys.path.append(os.path.abspath(os.path.dirname(__file__))+'/HXMT_GPS_py3Tools')
import multiprocessing as mp
from Lctools_201903 import FindLC_basefile
import time

def deal_newdata(lists):
    print('Generate LC:::',lists)
    ###Run_pool(fucgenlc,lists)
    with open(genlc_logfile,'a')as f:
        for i in lists:
            tmptm = str(time.strftime('%Y.%m.%d-%H:%M',time.localtime(time.time())))
            print(i,'Processed!',tmptm,file=f)

    with open(mkdir_path+'/gps_mission/all_mission.sh','a+')as f:
        for i in lists:
            print(
                'hep_sub /sharefs/hbkg/user/luoqi/GRB/work/ihep4/HXMT_GPS_Program/task_genlc.sh -argu {:s} -g hxmt -mem 8000 -o {:s} -e {:s}'.format(
                    i, scan_dir + '/joboutput', scan_dir + '/joboutput'), file=f)

    norm_lc = []
    try:
        with open(normlc_logtmp,'r')as f:
            for i in f:
                norm_lc.append(i[:13])
        print('rm %s'%normlc_logtmp)
        ###os.system('rm %s'%normlc_logtmp)
        with open(normlc_logfile,'a')as f:
            for i in norm_lc:
                tmptm = str(time.strftime('%Y.%m.%d-%H:%M',time.localtime(time.time())))
                print(i,' Normal-Generated!',tmptm,file=f)
    except:
        pass
    print(norm_lc)
    ###norm_lclist = [i[5:13] for i in norm_lc]
    if len(norm_lc)>0:
        print('Begin to Fit>>>')
        ###Run_pool(fuc_normfit,norm_lclist)


def run():
    folder_now = os.listdir(data1L_path)
    genlc_loglist = []
    with open(genlc_logfile,'a+')as f:
        for i in f:genlc_loglist.append(i[:11])
        #print(i[:11],i)
    print('ok')

    mins = set(folder_now)-set(genlc_loglist)
    folder_acess = []
    for i in mins:
        print(i)
        try:
            temp = FindLC_basefile(i+'01',data1L_path+i+'/',INST='HE',EHKFILE=False)
            folder_acess.append(i+'01')
            print('New folder accessible',i)
        except:
            print('New folder without enough FILES')
            pass
        #folder_acess.sort(reverse=True)
        #print folder_acess
    if len(folder_acess)>0:
        folder_acess.sort(reverse=True)
        print('New Obervation Accessible!',folder_acess)
        deal_newdata(folder_acess)
    else:
        print('There is no New Obervation!')
        time.sleep(1800)



if __name__ == "__main__":
    data1L_path = "/hxmt/work/HXMT-DATA/1L/A01/P0101295/"  # "/hxmt/work/HXMT-DATA/1L/A02/P0201013/"
    genlc_logfile = 'Genlc_process.log'
    normlc_logtmp = 'Genlc_tmp.log'
    normlc_logfile = 'Genlc_succes.log'
    srpath = '/sharefs/hbkg/users/luoqi/GRB/work/ihep4/data/HE_GPS_SML/'  # '/sharefs/hbkg/users/luoqi/GRB/work/ihep4/data/HE_GPS_SML/'


    dir_strlist = []
    with open('dir_config.txt', 'r', encoding='utf-8') as f:
        for line in f.readlines():
            line = line.strip()
            dir_strlist.append(line)
    print(dir_strlist)
    mkdir_path = dir_strlist[0]
    scan_dir = dir_strlist[1]
    #while True:
    run()

