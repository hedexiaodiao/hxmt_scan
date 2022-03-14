import numpy as np
import sys,os

if __name__ == "__main__":

    listall = []
    with open('A01_crab.txt','r')as f:
        for i in f:
            listall.append(i[:13])
    print(len(listall),listall)

    tot = len(listall)
    group_volum = 8
    group_all = int(tot/group_volum)
    group_remain = tot%group_volum
    print(group_all,group_remain,group_all*group_volum+group_remain)

    for i in range(group_all+1):
        if i<group_all:
            for j in range(group_volum):
                task_str = 'sh task_genlc.sh %s 7_20\n'%(listall[8*i+j])
                with open('task_A01_g%s.txt'% i, 'a+') as f:
                    print(task_str,file=f)
                calib_str = 'sh sub_gps_genlc.sh %s 7_20\n'%(listall[8*i+j])
                with open('calib_A01_g%s.txt'% i, 'a+') as f:
                    print(calib_str,file=f)
        if i==group_all:
            for j in range(group_remain):
                task_str = 'sh task_genlc.sh %s 7_20\n' % (listall[8 * i + j])
                with open('task_A01_g%s.txt' % i, 'a+') as f:
                    print(task_str,file=f)
                calib_str = 'sh sub_gps_genlc.sh %s 7_20\n' % (listall[8 * i + j])
                with open('calib_A01_g%s.txt' % i, 'a+') as f:
                    print(calib_str,file=f)