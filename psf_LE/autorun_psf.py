import numpy as np
import time
import sys,os
import glob
from xml.etree.ElementTree import ElementTree,Element
'''
本脚本将依次实现以下功能：
1.获取当前目录下的data目录列表，默认三层结构后是所需数据，即./data/*/*/box0.fits；
2.读取所需数据后，复制config.xml文件至1中所寻目录，将config里对应文件名更改为1中文件名；
3.写下当次将要执行的标定作业执行命令，至sh文件，可一键提交执行。
'''

#----------------get fits name---------------------#
def get_fitsname(fits_dir):
    att_list = glob.glob(r'%s/*att.fits' % (fits_dir))
    lc_list =  glob.glob(r'%s/*box*.fits' % (fits_dir))
    print(att_list,lc_list)
    return att_list,lc_list

#---------------make config------------------------#
def gen_config(att_list,lc_list,fits_dir,xmlpath='./config_he.xml',instru='LE'):
    exec_name = './command_'+instru.lower()+'_'+time.strftime("%y%m%d")+'.sh'

    tree = ElementTree()
    tree.parse(xmlpath)
    pathnodes = tree.findall("PATH/outpath")
    outnodes = tree.findall("Outfile_Name/outfilename")
    filenodes = tree.findall("Infile_Name/infilename/infilenamelist")
    instrunodes = tree.findall("Inst_Info/Instrument")
    outfile = fits_dir + "/config_%s.xml"%(instru.lower())
    print(outfile)
    filenodes[0].text = "\n  " + lc_list[0]+" \n"
    filenodes[1].text = "\n  " + lc_list[1]+" \n"
    filenodes[2].text = "\n  " + lc_list[2]+" \n"
    filenodes[3].text = '\n  ' + att_list[0] + ' \n'
    pathnodes[0].text = '\n  ' + fits_dir + '  \n'
    instrunodes[0].text = '\n ' +instru + ' \n'
    outnodes[0].text = '\n '+'/GAL_he_small \n'

    tree.write(outfile, encoding="utf-8", xml_declaration=True)
    with open(exec_name,'a+')as f:
        try:
            print('sh ./sub_task_psf.sh %s'%(outfile),file=f)
        except:
            print('cannot write down command!!!')


#-------------------count dir----------------------#
def count_dir(data_dir):
    tree1 = os.listdir(data_dir)
    i=0
    dir_list = []
    #print(tree1)
    for element in tree1:
        path = os.path.join(data_dir,element)
        #print(path)
        if os.path.isdir(path):
            i=i+1
            dir_list.append(path)
    print(dir_list,i)
    return dir_list,i


data_dir = './data'
if len(sys.argv)<2:
    print("Use the default Instrument LE!")
    sys.exit(1)
print(sys.argv[0],sys.argv[1])
instru = sys.argv[1]
tree1_dirlist, tree1_tot = count_dir(data_dir)
for i in range(tree1_tot):
    tree2_dirlist, tree2_tot = count_dir(tree1_dirlist[i])
    for j in range(tree2_tot):
        fits_dir = tree2_dirlist[j]
        att_list, lc_list = get_fitsname(fits_dir)
        gen_config(att_list,lc_list,fits_dir,instru=instru)