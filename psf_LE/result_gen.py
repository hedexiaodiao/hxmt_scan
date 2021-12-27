import numpy as np
import pandas as pd
import time,sys

'''
本脚本将依次实现以下功能：
读取提供的Para_log.txt文件，将其处理成二次型PSF模型所需输入参数
所需参数：
logtxt, 对应Para_log.txt的文件名
instrument, 用以命名生成文件所需, 代表所用仪器, 手动修改文件名也可
datayear, 用以命名生成文件所需, 代表所用数据年份, 手动修改文件名也可
txtname, 直接定义生成文件的全名, 可指定绝对路径
'''

def result_gen(logtxt='Para_log.txt',instrument='HE',datayear=2020,txtname=None):
    log_pd = pd.read_table(logtxt, sep=' ')
    log_np = log_pd.values
    result_box = np.array([])
    for i in range(1+0,1+3):
        p = 68*i - 12
        q = 68*i - 3
        result_box = np.append(result_box,log_np[p:q,2])
    arr_txt = np.ones((27,2))
    arr_txt[:,0] = result_box
    if txtname==None:
        txtname = 'Parabolic_'+instrument.upper()+'_'+str(datayear)+'_'+time.strftime("%y%m")+'.txt'
    np.savetxt(txtname,arr_txt)

if __name__ == "__main__":
    if len(sys.argv)<2:
        print("Use the default Para_log.txt!")
        print("Use the default Instrument HE!")
        print("Use the default data year 2020!")
        logtxt = 'Para_log.txt'
        instrument = 'HE'
        datayear = 2020
        txtname = None
    elif len(sys.argv)==2:
        print(sys.argv[0],sys.argv[1])
        print("Use the default Instrument HE!")
        print("Use the default data year 2020!")
        logtxt = sys.argv[1]
        instrument = 'HE'
        datayear = 2020
        txtname = None
    elif len(sys.argv)==3:
        print(sys.argv[0], sys.argv[1],sys.argv[2])
        print("Use the default data year 2020!")
        logtxt = sys.argv[1]
        instrument = sys.argv[2]
        datayear = 2020
        txtname = None
    elif len(sys.argv)==4:
        print(sys.argv[0], sys.argv[1],sys.argv[2],sys.argv[3])
        logtxt = sys.argv[1]
        instrument = sys.argv[2]
        datayear = sys.argv[3]
        txtname = None
    else:
        logtxt = sys.argv[1]
        instrument = sys.argv[2]
        datayear = sys.argv[3]
        txtname = sys.argv[4]

    result_gen(logtxt=logtxt, instrument=instrument, datayear=datayear, txtname=txtname)

