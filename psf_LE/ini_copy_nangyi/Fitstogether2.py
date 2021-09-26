import numpy as np
from astropy.io import fits as pf
import sys
import os
from xml.etree.ElementTree import ElementTree, Element
sys.path.append('/hxmt/work/USERS/nangyi/tools')
import Findfile
run = Findfile.findfile()

def fitstogether2(openfile, outpath):
    hd1 = pf.open(openfile[0])
    hd2 = pf.open(openfile[1])
    prihdr = pf.Header()
    prihdu = pf.PrimaryHDU(header = prihdr)
    hdulist = pf.HDUList([prihdu])
    print openfile, len(hd1) - 1
    for i in range(len(hd1) - 1):
        nr1 = hd1[i + 1].data.shape[0]
        nr2 = hd2[i + 1].data.shape[0]
        nr = nr1 + nr2
        hdu = pf.BinTableHDU.from_columns(hd1[i + 1].columns, nrows = nr,header=prihdr)
        for colname in hd1[i + 1].columns.names:
            hdu.data[colname][nr1:] = hd2[i + 1].data[colname]
        
        hdu.name = hd1[i + 1].name
        hdulist.append(hdu)
    
    hdulist.writeto(outpath, clobber = True)
    hd1.close()
    hd2.close()
    del hd1
    del hd2


def fitstogetherMulti(infilelist, outname):
    if len(infilelist) < 2:
        print 'lenth of input list should > 2!'
        raise Exception('lenth of input list should > 2!')
    fitstogether2([infilelist[0],infilelist[1]], outname)
    for i in infilelist[2:]:
        fitstogether2([outname,i], outname)
    
if __name__ == '__main__':
    if len(sys.argv) > 2:
        key1 = sys.argv[1]
        key2 = sys.argv[2]
    print key1, key2
    att1 = run(key1[:-2] + '_Att_FFFFFF_V1_1RP', '/hxmt/work/HXMT-DATA/1R/PIPROD/A01/P0101294/')
    att2 = run(key2[:-2] + '_Att_FFFFFF_V1_1RP', '/hxmt/work/HXMT-DATA/1R/PIPROD/A01/P0101294/')
    filestr = '/hxmt/work/HXMT_scan_data/HE/nyRt/' + key1[5:] + '_' + key2[-4:]
    
    try:
        os.makedirs(filestr, mode = 484)
    except OSError:
        pass
    filepath = [att1[0],att2[0]]
    print filepath
    print filestr
    outpathatt = filestr + '/' + key1[5:] + '_' + key2[-4:] + 'Att.fits'
    fitstogether2(filepath, outpathatt)
    datapath = '/hxmt/work/HXMT_scan_data/HE/data/'
    for i in [0,1,5]:
        filepath = [datapath + 'bkg/P0101294' + key1[-5:] + 'he_lc%s_bkg.fits' % i,datapath + 'bkg/P0101294' + key2[-5:] + 'he_lc%s_bkg.fits' % i]
        outpath = filestr + '/' + key1[5:] + '_' + key2[-4:] + 'he_lc%s_bkg.fits' % i
        fitstogether2(filepath, outpath)
    
    in_path = '/hxmt/work/HXMT_scan_data/HE/nangyi/config_he_small.xml'
    tree = ElementTree()
    tree.parse(in_path)
    pathnodes = tree.findall('PATH/outpath')
    outnodes = tree.findall('Outfile_Name/outfilename')
    filenodes = tree.findall('Infile_Name/infilename/infilenamelist')
    outfile = filestr + '/config_he.xml'
    print outfile
    filenodes[0].text = '\n  ' + outpath[:-15] + 'he_lc0_bkg.fits \n'
    filenodes[1].text = '\n  ' + outpath[:-15] + 'he_lc1_bkg.fits \n'
    filenodes[2].text = '\n  ' + outpath[:-15] + 'he_lc5_bkg.fits \n'
    filenodes[3].text = '\n  ' + outpathatt + ' \n'
    pathnodes[0].text = '\n  ' + filestr + '  \n'
    tree.write(outfile, encoding = 'utf-8', xml_declaration = True)
    os.system('cp loadmod_he2.py %s' % (filestr + '/'))
