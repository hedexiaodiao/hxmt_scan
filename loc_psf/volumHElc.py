#!/home/hxmt/nangyi/anaconda2/bin/python
#from multiprocessing import Process
import multiprocessing as mp
import sys
import os

listall = []
with open('crabpoint.txt','r')as f:
    for i in f:
        listall.append(i[:13])

print listall
def fuc(path):
    #cmd = 'python genhe.py %s note.xml'%path
    #cmd= 'python genme.py %s'%path
    cmd= 'python metest.py %s'%path
    print cmd
    temp = os.system(cmd)
    return [path,temp]

pool = mp.Pool(processes = 10)
res = pool.map(fuc, listall)


pool.close()
pool.join()
#for i in res:
#    with open('log_HELC.txt','a')as f:
#        f.write(i[0]+'\t%s\n'%i[1])


'''
if __name__ == '__main__':
    plog = []
    for i in listall:
        p = Process(target = fuc, args = (i,))
        p.start()
        plog.append(p)
    for p in plog:
        p.join()
'''
