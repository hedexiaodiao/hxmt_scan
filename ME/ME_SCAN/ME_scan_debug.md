##Notes
比较标准的HXMTSoft环境：
/home/hxmt/hxmtsoft2/hxmtsoftv2.04/hxmtsoft_v2.04.sh
A04目录结构与之前完全不一样，没有找到EHK文件
A03, /sharefs/hbkg/data/SCAN/ME/Org/GTI/P030124061101_gtiv2.fits是0k <---1L的EHK文件里 SUN_ANG, MOON_ANG从某时开始全是0

##A03
```buildoutcfg
megrade : HXMT ME task, megrade is running
megrade : PILParSet Warning: parameter 'clobber' set to yes!
megrade : PILParSet Warning: the file /sharefs/hbkg/data/SCAN/ME/Org/P030124061101/P030124061101_me_grade.fits will be overwritten!
[####################################################################################################][100%]
mepical : HXMT ME task, mepical is running successfully!
mepical : ##############################################
P030124061201  : pi end
megrade : ##############################################
megrade : HXMT ME task, megrade is running
megrade : PILParSet Warning: parameter 'clobber' set to yes!
megrade : PILParSet Warning: the file /sharefs/hbkg/data/SCAN/ME/Org/P030124061201/P030124061201_me_grade.fits will be overwritten!
[####################################################################################################][100%]
mepical : HXMT ME task, mepical is running successfully!
mepical : ##############################################
P030124043001  : pi end
megrade : ##############################################
megrade : HXMT ME task, megrade is running
megrade : PILParSet Warning: parameter 'clobber' set to yes!
megrade : PILParSet Warning: the file /sharefs/hbkg/data/SCAN/ME/Org/P030124043001/P030124043001_me_grade.fits will be overwritten!
[####################################################################################################][100%]
megrade : PILParSet Warning: parameter 'clobber' set to yes!
megrade : PILParSet Warning: the file /sharefs/hbkg/data/SCAN/ME/Org/P030124013801/P030124013801_me_dead.fits will be overwritten!
megrade : HXMT ME task, megrade is running successfully!
megrade : ##############################################
P030124013801  : grade end
rm: cannot remove ‘/sharefs/hbkg/data/SCAN/ME/Org/GTI/P030124013801_me_gti.fits’: No such file or directory
megtigen : ##############################################                                            ][39%]
megtigen : HXMT ME task, megtigen is running
megtigen : HXMT ME task, megtigen is running successfully!
megtigen : ##############################################
P030124013801  : o-gti & ehk end
Process Process-12:
Traceback (most recent call last):
  File "/hxmt/soft/Develop/anaconda2/lib/python2.7/multiprocessing/process.py", line 267, in _bootstrap
    self.run()
  File "/hxmt/soft/Develop/anaconda2/lib/python2.7/multiprocessing/process.py", line 114, in run
    self._target(*self._args, **self._kwargs)
  File "multi_time_merun.py", line 22, in fuc
    timing_run.merun_v2(fpath,Wpath,str(path))
  File "/sharefs/hbkg/user/luoqi/HXMT_SCAN/ME/ME_SCAN/genlc/timing_run.py", line 191, in merun_v2
    regti_v2.regti(meevtfile, Wpath, ObsID)###calc total time
  File "/sharefs/hbkg/user/luoqi/HXMT_SCAN/ME/ME_SCAN/genlc/regti_v2.py", line 25, in regti
    evt_data = evt_all[1].data
  File "/hxmt/soft/Develop/anaconda2/lib/python2.7/site-packages/astropy/utils/decorators.py", line 739, in __get__
    val = self.fget(obj)
  File "/hxmt/soft/Develop/anaconda2/lib/python2.7/site-packages/astropy/io/fits/hdu/table.py", line 350, in data
    data = self._get_tbdata()
  File "/hxmt/soft/Develop/anaconda2/lib/python2.7/site-packages/astropy/io/fits/hdu/table.py", line 155, in _get_tbdata
    columns = self.columns
  File "/hxmt/soft/Develop/anaconda2/lib/python2.7/site-packages/astropy/utils/decorators.py", line 739, in __get__
    val = self.fget(obj)
  File "/hxmt/soft/Develop/anaconda2/lib/python2.7/site-packages/astropy/io/fits/hdu/table.py", line 346, in columns
    return self._columns_type(self)
  File "/hxmt/soft/Develop/anaconda2/lib/python2.7/site-packages/astropy/io/fits/column.py", line 1195, in __init__
    self._init_from_table(input)
  File "/hxmt/soft/Develop/anaconda2/lib/python2.7/site-packages/astropy/io/fits/column.py", line 1259, in _init_from_table
    for keyword, value in iteritems(hdr):
  File "/hxmt/soft/Develop/anaconda2/lib/python2.7/site-packages/astropy/io/fits/header.py", line 940, in iteritems
    yield (card.keyword, card.value)
  File "/hxmt/soft/Develop/anaconda2/lib/python2.7/site-packages/astropy/io/fits/card.py", line 285, in value
    self._value = self._parse_value()
  File "/hxmt/soft/Develop/anaconda2/lib/python2.7/site-packages/astropy/io/fits/card.py", line 744, in _parse_value
    ".verify('fix').".format(self.keyword))
VerifyError: Unparsable card (OBS_ID), fix it first with .verify('fix').
[####################################################################################################][100%]
megrade : PILParSet Warning: parameter 'clobber' set to yes!
megrade : PILParSet Warning: the file /sharefs/hbkg/data/SCAN/ME/Org/P030124061101/P030124061101_me_dead.fits will be overwritten!
megrade : HXMT ME task, megrade is running successfully!
megrade : ##############################################
P030124061101  : grade end
rm: cannot remove ‘/sharefs/hbkg/data/SCAN/ME/Org/GTI/P030124061101_me_gti.fits’: No such file or directory
megtigen : ########################################################################################  ][98%]
megtigen : HXMT ME task, megtigen is running
megtigen : GtiCalc warning: no GTI extension in GTI List!#######################################     ][95%]
megtigen : GtiCalc: Can't create GTI File!
megtigen : HXMT ME task, megtigen is running unsuccessfully!
megtigen : ##############################################
P030124061101  : o-gti & ehk end
Process Process-14:
Traceback (most recent call last):
  File "/hxmt/soft/Develop/anaconda2/lib/python2.7/multiprocessing/process.py", line 267, in _bootstrap
    self.run()
  File "/hxmt/soft/Develop/anaconda2/lib/python2.7/multiprocessing/process.py", line 114, in run
    self._target(*self._args, **self._kwargs)
  File "multi_time_merun.py", line 22, in fuc
    timing_run.merun_v2(fpath,Wpath,str(path))
  File "/sharefs/hbkg/user/luoqi/HXMT_SCAN/ME/ME_SCAN/genlc/timing_run.py", line 191, in merun_v2
    regti_v2.regti(meevtfile, Wpath, ObsID)###calc total time
  File "/sharefs/hbkg/user/luoqi/HXMT_SCAN/ME/ME_SCAN/genlc/regti_v2.py", line 22, in regti
    soft_gti = pf.open(Wgtifile)
  File "/hxmt/soft/Develop/anaconda2/lib/python2.7/site-packages/astropy/io/fits/hdu/hdulist.py", line 148, in fitsopen
    lazy_load_hdus, **kwargs)
  File "/hxmt/soft/Develop/anaconda2/lib/python2.7/site-packages/astropy/io/fits/hdu/hdulist.py", line 402, in fromfile
    lazy_load_hdus=lazy_load_hdus, **kwargs)
  File "/hxmt/soft/Develop/anaconda2/lib/python2.7/site-packages/astropy/io/fits/hdu/hdulist.py", line 1040, in _readfrom
    raise IOError('Empty or corrupt FITS file')
IOError: Empty or corrupt FITS file
[####################################################################################################][100%]
megrade : PILParSet Warning: parameter 'clobber' set to yes!
megrade : PILParSet Warning: the file /sharefs/hbkg/data/SCAN/ME/Org/P030124061201/P030124061201_me_dead.fits will be overwritten!
megrade : HXMT ME task, megrade is running successfully!
megrade : ##############################################
P030124061201  : grade end
rm: cannot remove ‘/sharefs/hbkg/data/SCAN/ME/Org/GTI/P030124061201_me_gti.fits’: No such file or directory
megtigen : ##############################################
megtigen : HXMT ME task, megtigen is running
megtigen : GtiCalc warning: no GTI extension in GTI List!
megtigen : GtiCalc: Can't create GTI File!
megtigen : HXMT ME task, megtigen is running unsuccessfully!
megtigen : ##############################################
P030124061201  : o-gti & ehk end
Process Process-15:
Traceback (most recent call last):
  File "/hxmt/soft/Develop/anaconda2/lib/python2.7/multiprocessing/process.py", line 267, in _bootstrap
    self.run()
  File "/hxmt/soft/Develop/anaconda2/lib/python2.7/multiprocessing/process.py", line 114, in run
    self._target(*self._args, **self._kwargs)
  File "multi_time_merun.py", line 22, in fuc
    timing_run.merun_v2(fpath,Wpath,str(path))
  File "/sharefs/hbkg/user/luoqi/HXMT_SCAN/ME/ME_SCAN/genlc/timing_run.py", line 191, in merun_v2
    regti_v2.regti(meevtfile, Wpath, ObsID)###calc total time
  File "/sharefs/hbkg/user/luoqi/HXMT_SCAN/ME/ME_SCAN/genlc/regti_v2.py", line 22, in regti
    soft_gti = pf.open(Wgtifile)
  File "/hxmt/soft/Develop/anaconda2/lib/python2.7/site-packages/astropy/io/fits/hdu/hdulist.py", line 148, in fitsopen
    lazy_load_hdus, **kwargs)
  File "/hxmt/soft/Develop/anaconda2/lib/python2.7/site-packages/astropy/io/fits/hdu/hdulist.py", line 402, in fromfile
    lazy_load_hdus=lazy_load_hdus, **kwargs)
  File "/hxmt/soft/Develop/anaconda2/lib/python2.7/site-packages/astropy/io/fits/hdu/hdulist.py", line 1040, in _readfrom
    raise IOError('Empty or corrupt FITS file')
IOError: Empty or corrupt FITS file
[####################################################################################################][100%]
megrade : PILParSet Warning: parameter 'clobber' set to yes!
megrade : PILParSet Warning: the file /sharefs/hbkg/data/SCAN/ME/Org/P030124043001/P030124043001_me_dead.fits will be overwritten!
megrade : HXMT ME task, megrade is running successfully!
megrade : ##############################################
P030124043001  : grade end
megtigen : ##############################################
megtigen : HXMT ME task, megtigen is running
megtigen : HXMT ME task, megtigen is running successfully!
megtigen : ##############################################
P030124043001  : o-gti & ehk end
[2.87328129e+08] [2.87339923e+08]
WARNING: AstropyDeprecationWarning: "clobber" was deprecated in version 2.0 and will be removed in a future version. Use argument "overwrite" instead. [astropy.utils.decorators]
P030124043001  : r-gti end
mescreen : ##############################################
mescreen : HXMT ME task, mescreen is running
mescreen : PILParSet Warning: parameter 'clobber' set to yes!
mescreen : PILParSet Warning: the file /sharefs/hbkg/data/SCAN/ME/Org/Screen/P030124043001_me_screen.fits will be overwritten!
mescreen : mescreen: User detector selection:
mescreen : Bad (Bad,hot,flickening) detector List:
mescreen :        ID:  1

```

## A04
```buildoutcfg
megrade : HXMT ME task, megrade is running successfully!
megrade : ##############################################
P040124006801  : grade end
Process Process-92:
Traceback (most recent call last):
  File "/hxmt/soft/Develop/anaconda2/lib/python2.7/multiprocessing/process.py", line 267, in _bootstrap
    self.run()
  File "/hxmt/soft/Develop/anaconda2/lib/python2.7/multiprocessing/process.py", line 114, in run
    self._target(*self._args, **self._kwargs)
  File "multi_time_merun.py", line 22, in fuc
    timing_run.merun_v2(fpath,Wpath,str(path))
  File "/sharefs/hbkg/user/luoqi/HXMT_SCAN/ME/ME_SCAN/genlc/timing_run.py", line 178, in merun_v2
    os.system("cp %s %s" % (ehkfile, Wehkfile))
UnboundLocalError: local variable 'ehkfile' referenced before assignment
[####################################################################################################][100%]
megrade : HXMT ME task, megrade is running successfully!
megrade : ##############################################
P040124006901  : grade end
Process Process-93:
Traceback (most recent call last):
  File "/hxmt/soft/Develop/anaconda2/lib/python2.7/multiprocessing/process.py", line 267, in _bootstrap
    self.run()
  File "/hxmt/soft/Develop/anaconda2/lib/python2.7/multiprocessing/process.py", line 114, in run
    self._target(*self._args, **self._kwargs)
  File "multi_time_merun.py", line 22, in fuc
    timing_run.merun_v2(fpath,Wpath,str(path))
  File "/sharefs/hbkg/user/luoqi/HXMT_SCAN/ME/ME_SCAN/genlc/timing_run.py", line 178, in merun_v2
    os.system("cp %s %s" % (ehkfile, Wehkfile))
UnboundLocalError: local variable 'ehkfile' referenced before assignment
[####################################################################################################][100%]
megrade : HXMT ME task, megrade is running successfully!
megrade : ##############################################
P040124007001  : grade end
Process Process-94:
Traceback (most recent call last):
  File "/hxmt/soft/Develop/anaconda2/lib/python2.7/multiprocessing/process.py", line 267, in _bootstrap
    self.run()
  File "/hxmt/soft/Develop/anaconda2/lib/python2.7/multiprocessing/process.py", line 114, in run
    self._target(*self._args, **self._kwargs)
  File "multi_time_merun.py", line 22, in fuc
    timing_run.merun_v2(fpath,Wpath,str(path))
  File "/sharefs/hbkg/user/luoqi/HXMT_SCAN/ME/ME_SCAN/genlc/timing_run.py", line 178, in merun_v2
    os.system("cp %s %s" % (ehkfile, Wehkfile))
UnboundLocalError: local variable 'ehkfile' referenced before assignment
Check END
```

##After rm MOON/SUN ANG
####errorCorrect 56
由bdet未正常生成导致
```buildoutcfg
mescreen evtfile=/sharefs/hbkg/data/SCAN/ME/Org/P030124065601/P030124065601_me_grade.fits gtifile=/sharefs/hbkg/data/SCAN/ME/Org/GTI/P030124065601_me_gti.fits outfile=/sharefs/hbkg/data/SCAN/ME/Org/Screen/P030124065601_me_screen.fits baddetfile=/sharefs/hbkg/data/SCAN/ME/Org/P030124065601/P030124065601_me_bdet.fits userdetid="0-53"
mescreen : ##############################################
mescreen : HXMT ME task, mescreen is running
mescreen : mescreen: User detector selection:
mescreen : mescreen: Error: Unable to get 'merangefile'! Always in $MEADAS/refdata.
mescreen : HXMT ME task, mescreen is running unsuccessfully!
mescreen : ##############################################
```
```buildoutcfg
melcgen : HXMT ME task, melcgen is running unsuccessfully!
melcgen : ##############################################
P030124065601  : o-lc end
P030124065601  1.Soft Cal: End##############################                                         ][59%]
Process Process-24:
Traceback (most recent call last):
  File "/hxmt/soft/Develop/anaconda2/lib/python2.7/multiprocessing/process.py", line 267, in _bootstrap
    self.run()
  File "/hxmt/soft/Develop/anaconda2/lib/python2.7/multiprocessing/process.py", line 114, in run
    self._target(*self._args, **self._kwargs)
  File "multi_time_merun.py", line 22, in fuc
    timing_run.merun_v2(fpath,Wpath,str(path))
  File "/sharefs/hbkg/user/luoqi/HXMT_SCAN/ME/ME_SCAN/genlc/timing_run.py", line 241, in merun_v2
    errorCorrect.errorCorrect('%s' % Wpath,'%s' % ObsID)
  File "/sharefs/hbkg/user/luoqi/HXMT_SCAN/ME/ME_SCAN/genlc/errorCorrect.py", line 56, in errorCorrect
    data1=pf.open(slcfile_list[num])#########################                                        ][60%]
IndexError: list index out of range
[####################################################################################################][100%]
melcgen : HXMT ME task, melcgen is running successfully!
melcgen : ##############################################
P030124061801  : o-lc end
P030124061801  1.Soft Cal: End
```

####gen_lc_me 48
仍是SUN_ANG, MOON_ANG导致的,修正后解决这部分
```buildoutcfg
[   0 4546 8930] pot:: 1
::::: 4380 4380
1
1
1
56.0 -29.0
WARNING: AstropyDeprecationWarning: "clobber" was deprecated in version 2.0 and will be removed in a future version. Use argument "overwrite" instead. [astropy.utils.decorators]
WARNING: AstropyDeprecationWarning: "clobber" was deprecated in version 2.0 and will be removed in a future version. Use argument "overwrite" instead. [astropy.utils.decorators]
me 2 done
P030124061801 2.BKG Cal: LC without BKG and BKG value.
Process Process-21:
Traceback (most recent call last):
  File "/hxmt/soft/Develop/anaconda2/lib/python2.7/multiprocessing/process.py", line 267, in _bootstrap
    self.run()
  File "/hxmt/soft/Develop/anaconda2/lib/python2.7/multiprocessing/process.py", line 114, in run
    self._target(*self._args, **self._kwargs)
  File "multi_time_merun.py", line 22, in fuc
    timing_run.merun_v2(fpath,Wpath,str(path))
  File "/sharefs/hbkg/user/luoqi/HXMT_SCAN/ME/ME_SCAN/genlc/timing_run.py", line 247, in merun_v2
    gen_lc_me.gen_lc_me(Wehkfile, Wme_gtifile, '%s' % (Wpath), "%s" % ObsID)
  File "/sharefs/hbkg/user/luoqi/HXMT_SCAN/ME/ME_SCAN/genlc/gen_lc_me.py", line 48, in gen_lc_me
    gtitstart = np.r_[gti_all[0],gti_start]
IndexError: index 0 is out of bounds for axis 0 with size 0
P040124000201
/hxmt/work/HXMT-DATA/1L/A03/P0301240 P040124000201 : No data
P040124000301
```

###After fix 3 errors in ME_SCAN
P030124063901、P030124064001
```buildoutcfg
P030124064001  : grade end
rm: cannot remove ‘/sharefs/hbkg/data/SCAN/ME/Org/GTI/P030124064001_me_gti.fits’: No such file or directory
megtigen tempfile=/hxmt/work/HXMT-DATA/1L/A03/P0301240/P0301240640/P030124064001-20210706-01-01/ME/HXMT_P030124064001_ME-TH_FFFFFF_V1_L1P.FITS ehkfile=/sharefs/hbkg/data/SCAN/ME/Org/EHK/P030124064001_ehk.fits outfile=/sharefs/hbkg/data/SCAN/ME/Org/GTI/P030124064001_gtiv2.fits defaultexpr=NONE expr="ELV>5&&COR>=8&&T_SAA>=200&&TN_SAA>=100&&SAA_FLAG==0&&ANG_DIST<=359&&(SAT_LAT<31||SAT_LAT>38)&&(SAT_LON>245||SAT_LON<228)&&(SAT_LAT>=-36.5&&SAT_LAT<=36.5)"
megtigen : ##############################################
megtigen : HXMT ME task, megtigen is running
megtigen : GtiCalc warning: no GTI extension in GTI List!
megtigen : GtiCalc: Can't create GTI File!
megtigen : HXMT ME task, megtigen is running unsuccessfully!
megtigen : ##############################################
P030124064001  : o-gti & ehk end
Process Process-22:
Traceback (most recent call last):
  File "/hxmt/soft/Develop/anaconda2/lib/python2.7/multiprocessing/process.py", line 267, in _bootstrap
    self.run()
  File "/hxmt/soft/Develop/anaconda2/lib/python2.7/multiprocessing/process.py", line 114, in run
    self._target(*self._args, **self._kwargs)
  File "multi_time_merun.py", line 22, in fuc
    timing_run.merun_v2(fpath,Wpath,str(path))
  File "/sharefs/hbkg/user/luoqi/HXMT_SCAN/ME/ME_SCAN/genlc/timing_run.py", line 197, in merun_v2
    regti_v2.regti(meevtfile, Wpath, ObsID)###calc total time
  File "/sharefs/hbkg/user/luoqi/HXMT_SCAN/ME/ME_SCAN/genlc/regti_v2.py", line 22, in regti
    soft_gti = pf.open(Wgtifile)
  File "/hxmt/soft/Develop/anaconda2/lib/python2.7/site-packages/astropy/io/fits/hdu/hdulist.py", line 148, in fitsopen
    lazy_load_hdus, **kwargs)
  File "/hxmt/soft/Develop/anaconda2/lib/python2.7/site-packages/astropy/io/fits/hdu/hdulist.py", line 402, in fromfile
    lazy_load_hdus=lazy_load_hdus, **kwargs)
  File "/hxmt/soft/Develop/anaconda2/lib/python2.7/site-packages/astropy/io/fits/hdu/hdulist.py", line 1040, in _readfrom
    raise IOError('Empty or corrupt FITS file')
IOError: Empty or corrupt FITS file
```

###new
```buildoutcfg
megtigen : ##############################################                                            ][43%]
megtigen : HXMT ME task, megtigen is running
megtigen : HXMT ME task, megtigen is running successfully!
megtigen : ##############################################
P030124013801  : o-gti & ehk end
Process Process-12:
Traceback (most recent call last):
  File "/hxmt/soft/Develop/anaconda2/lib/python2.7/multiprocessing/process.py", line 267, in _bootstrap
    self.run()
  File "/hxmt/soft/Develop/anaconda2/lib/python2.7/multiprocessing/process.py", line 114, in run
    self._target(*self._args, **self._kwargs)
  File "multi_time_merun.py", line 22, in fuc
    timing_run.merun_v2(fpath,Wpath,str(path))
  File "/sharefs/hbkg/user/luoqi/HXMT_SCAN/ME/ME_SCAN/genlc/timing_run.py", line 197, in merun_v2
    regti_v2.regti(meevtfile, Wpath, ObsID)###calc total time
  File "/sharefs/hbkg/user/luoqi/HXMT_SCAN/ME/ME_SCAN/genlc/regti_v2.py", line 25, in regti
    evt_data = evt_all[1].data
  File "/hxmt/soft/Develop/anaconda2/lib/python2.7/site-packages/astropy/utils/decorators.py", line 739, in __get__
    val = self.fget(obj)
  File "/hxmt/soft/Develop/anaconda2/lib/python2.7/site-packages/astropy/io/fits/hdu/table.py", line 350, in data
    data = self._get_tbdata()
  File "/hxmt/soft/Develop/anaconda2/lib/python2.7/site-packages/astropy/io/fits/hdu/table.py", line 155, in _get_tbdata
    columns = self.columns
  File "/hxmt/soft/Develop/anaconda2/lib/python2.7/site-packages/astropy/utils/decorators.py", line 739, in __get__
    val = self.fget(obj)
  File "/hxmt/soft/Develop/anaconda2/lib/python2.7/site-packages/astropy/io/fits/hdu/table.py", line 346, in columns
    return self._columns_type(self)
  File "/hxmt/soft/Develop/anaconda2/lib/python2.7/site-packages/astropy/io/fits/column.py", line 1195, in __init__
    self._init_from_table(input)
  File "/hxmt/soft/Develop/anaconda2/lib/python2.7/site-packages/astropy/io/fits/column.py", line 1259, in _init_from_table
    for keyword, value in iteritems(hdr):
  File "/hxmt/soft/Develop/anaconda2/lib/python2.7/site-packages/astropy/io/fits/header.py", line 940, in iteritems
    yield (card.keyword, card.value)
  File "/hxmt/soft/Develop/anaconda2/lib/python2.7/site-packages/astropy/io/fits/card.py", line 285, in value
    self._value = self._parse_value()
  File "/hxmt/soft/Develop/anaconda2/lib/python2.7/site-packages/astropy/io/fits/card.py", line 744, in _parse_value
    ".verify('fix').".format(self.keyword))
VerifyError: Unparsable card (OBS_ID), fix it first with .verify('fix').
```

```buildoutcfg
  File "/hxmt/soft/Develop/anaconda2/lib/python2.7/site-packages/astropy/io/fits/column.py", line 1259, in _init_from_table
    for keyword, value in iteritems(hdr):
  File "/hxmt/soft/Develop/anaconda2/lib/python2.7/site-packages/astropy/io/fits/header.py", line 940, in iteritems
    yield (card.keyword, card.value)
  File "/hxmt/soft/Develop/anaconda2/lib/python2.7/site-packages/astropy/io/fits/card.py", line 285, in value
    self._value = self._parse_value()
  File "/hxmt/soft/Develop/anaconda2/lib/python2.7/site-packages/astropy/io/fits/card.py", line 744, in _parse_value
    ".verify('fix').".format(self.keyword))
VerifyError: Unparsable card (OBS_ID), fix it first with .verify('fix').
^CTraceback (most recent call last):########################################                         ][75%]
  File "multi_time_merun.py", line 174, in <module>
    main(0,0)
  File "multi_time_merun.py", line 137, in main
P030124064001  : grade end
P030124063901  : grade end
P030124043001  : grade end
    doSth(list5,fpath,lcpath)
  File "multi_time_merun.py", line 38, in doSth
    p.join()
  File "/hxmt/soft/Develop/anaconda2/lib/python2.7/multiprocessing/process.py", line 148, in join
    res = self._popen.wait(timeout)
  File "/hxmt/soft/Develop/anaconda2/lib/python2.7/multiprocessing/forking.py", line 154, in wait
    return self.poll(0)
  File "/hxmt/soft/Develop/anaconda2/lib/python2.7/multiprocessing/forking.py", line 135, in poll
    pid, sts = os.waitpid(self.pid, flag)
```

```buildoutcfg
megtigen : ##############################################
megtigen : GtiCalc warning: no GTI extension in GTI List!
megtigen : GtiCalc: Can't create GTI File!
megtigen : HXMT ME task, megtigen is running unsuccessfully!
megtigen : ##############################################
P030124063901  : o-gti & ehk end
P030124043001  : o-gti & ehk end
Process Process-14:
Traceback (most recent call last):
  File "/hxmt/soft/Develop/anaconda2/lib/python2.7/multiprocessing/process.py", line 267, in _bootstrap
    self.run()
  File "/hxmt/soft/Develop/anaconda2/lib/python2.7/multiprocessing/process.py", line 114, in run
    self._target(*self._args, **self._kwargs)
  File "multi_time_merun.py", line 22, in fuc
    timing_run.merun_v2(fpath,Wpath,str(path))
  File "/sharefs/hbkg/user/luoqi/HXMT_SCAN/ME/ME_SCAN/genlc/timing_run.py", line 197, in merun_v2
    regti_v2.regti(meevtfile, Wpath, ObsID)###calc total time
  File "/sharefs/hbkg/user/luoqi/HXMT_SCAN/ME/ME_SCAN/genlc/regti_v2.py", line 22, in regti
    soft_gti = pf.open(Wgtifile)
  File "/hxmt/soft/Develop/anaconda2/lib/python2.7/site-packages/astropy/io/fits/hdu/hdulist.py", line 148, in fitsopen
    lazy_load_hdus, **kwargs)
  File "/hxmt/soft/Develop/anaconda2/lib/python2.7/site-packages/astropy/io/fits/hdu/hdulist.py", line 402, in fromfile
    lazy_load_hdus=lazy_load_hdus, **kwargs)
  File "/hxmt/soft/Develop/anaconda2/lib/python2.7/site-packages/astropy/io/fits/hdu/hdulist.py", line 1040, in _readfrom
    raise IOError('Empty or corrupt FITS file')
IOError: Empty or corrupt FITS file
```

```buildoutcfg
P030124043001 
Process Process-13:
Traceback (most recent call last):
  File "/hxmt/soft/Develop/anaconda2/lib/python2.7/multiprocessing/process.py", line 267, in _bootstrap
    self.run()
  File "/hxmt/soft/Develop/anaconda2/lib/python2.7/multiprocessing/process.py", line 114, in run
    self._target(*self._args, **self._kwargs)
  File "multi_time_merun.py", line 22, in fuc
    timing_run.merun_v2(fpath,Wpath,str(path))
  File "/sharefs/hbkg/user/luoqi/HXMT_SCAN/ME/ME_SCAN/genlc/timing_run.py", line 244, in merun_v2
    nbkghesmtpoly2.poly_bkg('%s' % NetWpath,'%s' % Wpath,'%s' % ObsID)
  File "/sharefs/hbkg/user/luoqi/HXMT_SCAN/ME/ME_SCAN/genlc/nbkghesmtpoly2.py", line 274, in poly_bkg
    xr2_small,yr2_small,bkg2_small,error2 = get_bkg(slcfile_list[2])
  File "/sharefs/hbkg/user/luoqi/HXMT_SCAN/ME/ME_SCAN/genlc/nbkghesmtpoly2.py", line 195, in get_bkg
    tp = frag.min()
  File "/hxmt/soft/Develop/anaconda2/lib/python2.7/site-packages/numpy/core/_methods.py", line 32, in _amin
    return umr_minimum(a, axis, None, out, keepdims, initial)
ValueError: zero-size array to reduction operation minimum which has no identity
```

### use 1N ehk
```buildoutcfg
P030124043001  : pi end
megrade : ##############################################
megrade : HXMT ME task, megrade is running
megrade : PILParSet Warning: parameter 'clobber' set to yes!
megrade : PILParSet Warning: the file /sharefs/hbkg/data/SCAN/ME/Org/P030124043001/P030124043001_me_grade.fits will be overwritten!
[####################################################################################################][100%]
megrade : PILParSet Warning: parameter 'clobber' set to yes!
megrade : PILParSet Warning: the file /sharefs/hbkg/data/SCAN/ME/Org/P030124013801/P030124013801_me_dead.fits will be overwritten!
megrade : HXMT ME task, megrade is running successfully!
megrade : ##############################################
P030124013801  : grade end
Process Process-12:
Traceback (most recent call last):
  File "/hxmt/soft/Develop/anaconda2/lib/python2.7/multiprocessing/process.py", line 267, in _bootstrap
    self.run()
  File "/hxmt/soft/Develop/anaconda2/lib/python2.7/multiprocessing/process.py", line 114, in run
    self._target(*self._args, **self._kwargs)
  File "multi_time_merun.py", line 22, in fuc
    timing_run.merun_v2(fpath,Wpath,str(path))
  File "/sharefs/hbkg/user/luoqi/HXMT_SCAN/ME/ME_SCAN/genlc/timing_run.py", line 183, in merun_v2
    os.system("cp %s %s" % (ehkfile, Wehkfile))
UnboundLocalError: local variable 'ehkfile' referenced before assignment
```
```buildoutcfg
megrade : PILParSet Warning: parameter 'clobber' set to yes!
megrade : PILParSet Warning: the file /sharefs/hbkg/data/SCAN/ME/Org/P030124064001/P030124064001_me_dead.fits will be overwritten!
megrade : HXMT ME task, megrade is running successfully!
megrade : ##############################################
P030124064001  : grade end
rm: cannot remove ‘/sharefs/hbkg/data/SCAN/ME/Org/GTI/P030124064001_me_gti.fits’: No such file or directory
megtigen tempfile=/hxmt/work/HXMT-DATA/1L/A03/P0301240/P0301240640/P030124064001-20210706-01-01/ME/HXMT_P030124064001_ME-TH_FFFFFF_V1_L1P.FITS ehkfile=/sharefs/hbkg/data/SCAN/ME/Org/EHK/P030124064001_ehk.fits outfile=/sharefs/hbkg/data/SCAN/ME/Org/GTI/P030124064001_gtiv2.fits defaultexpr=NONE expr="ELV>5&&COR>=8&&T_SAA>=200&&TN_SAA>=100&&SAA_FLAG==0&&SUN_ANG>=10&&MOON_ANG>=5&&ANG_DIST<=359&&(SAT_LAT<31||SAT_LAT>38)&&(SAT_LON>245||SAT_LON<228)&&(SAT_LAT>=-36.5&&SAT_LAT<=36.5)"
megtigen : ##############################################
megtigen : HXMT ME task, megtigen is running
megtigen : HXMT ME task, megtigen is running successfully!############                               ][69%]
megtigen : ##############################################
P030124064001  : o-gti & ehk end
[3.00238379e+08] [3.00250048e+08]##########################################                          ][74%]
WARNING: AstropyDeprecationWarning: "clobber" was deprecated in version 2.0 and will be removed in a future version. Use argument "overwrite" instead. [astropy.utils.decorators]
P030124064001  : r-gti end
q list size:  (71060,)
291###########################################################################                       ][77%]
srcs::::::: [[147.33446  251.44875  -45.61111 ]
 [462.438    256.43542  -36.423058]
 [ 28.75739  262.98917  -33.83472 ]
 [ 96.02605  257.225    -44.102   ]]
gti sum:  0.0
Process Process-14:
Traceback (most recent call last):
  File "/hxmt/soft/Develop/anaconda2/lib/python2.7/multiprocessing/process.py", line 267, in _bootstrap
    self.run()
  File "/hxmt/soft/Develop/anaconda2/lib/python2.7/multiprocessing/process.py", line 114, in run
    self._target(*self._args, **self._kwargs)
  File "multi_time_merun.py", line 22, in fuc
    timing_run.merun_v2(fpath,Wpath,str(path))
  File "/sharefs/hbkg/user/luoqi/HXMT_SCAN/ME/ME_SCAN/genlc/timing_run.py", line 203, in merun_v2
    mescreen_newbd.mescreen_newbd('%s' % Wpath,'%s' % ObsID)###
  File "/sharefs/hbkg/user/luoqi/HXMT_SCAN/ME/ME_SCAN/genlc/mescreen_newbd.py", line 281, in mescreen_newbd
    new_str,new_stp=bkg_gti(Wpath, ObsID, ft)
TypeError: 'NoneType' object is not iterable
```
P030124043001
```buildoutcfg
Process Process-13:
Traceback (most recent call last):
  File "/hxmt/soft/Develop/anaconda2/lib/python2.7/multiprocessing/process.py", line 267, in _bootstrap
    self.run()
  File "/hxmt/soft/Develop/anaconda2/lib/python2.7/multiprocessing/process.py", line 114, in run
    self._target(*self._args, **self._kwargs)
  File "multi_time_merun.py", line 22, in fuc
    timing_run.merun_v2(fpath,Wpath,str(path))
  File "/sharefs/hbkg/user/luoqi/HXMT_SCAN/ME/ME_SCAN/genlc/timing_run.py", line 246, in merun_v2
    nbkghesmtpoly2.poly_bkg('%s' % NetWpath,'%s' % Wpath,'%s' % ObsID)
  File "/sharefs/hbkg/user/luoqi/HXMT_SCAN/ME/ME_SCAN/genlc/nbkghesmtpoly2.py", line 274, in poly_bkg
    xr2_small,yr2_small,bkg2_small,error2 = get_bkg(slcfile_list[2])
  File "/sharefs/hbkg/user/luoqi/HXMT_SCAN/ME/ME_SCAN/genlc/nbkghesmtpoly2.py", line 195, in get_bkg
    tp = frag.min()
  File "/hxmt/soft/Develop/anaconda2/lib/python2.7/site-packages/numpy/core/_methods.py", line 32, in _amin
    return umr_minimum(a, axis, None, out, keepdims, initial)
ValueError: zero-size array to reduction operation minimum which has no identity
```
multi-run bug
```buildoutcfg
P040124000201
P040124000301
/hxmt/work/HXMT-DATA/1L/A03/P0301240 P040124000201 : No data
P040124000401
P040124000501
P040124000601
/hxmt/work/HXMT-DATA/1L/A03/P0301240 P040124000301 : No data
/hxmt/work/HXMT-DATA/1L/A03/P0301240 P040124000401 : No data
/hxmt/work/HXMT-DATA/1L/A03/P0301240 P040124000501 : No data
/hxmt/work/HXMT-DATA/1L/A03/P0301240 P040124000601 : No data
P040124000701
P040124000801
P040124000901
P040124001001
/hxmt/work/HXMT-DATA/1L/A03/P0301240 P040124000701 : No data
P040124001101
/hxmt/work/HXMT-DATA/1L/A03/P0301240 P040124000801 : No data
/hxmt/work/HXMT-DATA/1L/A03/P0301240 P040124000901 : No data
/hxmt/work/HXMT-DATA/1L/A03/P0301240 P040124001001 : No data
/hxmt/work/HXMT-DATA/1L/A03/P0301240 P040124001101 : No data
P040124001201
P040124001301
P040124001401
P040124001501
/hxmt/work/HXMT-DATA/1L/A03/P0301240 P040124001201 : No data
P040124001601
/hxmt/work/HXMT-DATA/1L/A03/P0301240 P040124001301 : No data
/hxmt/work/HXMT-DATA/1L/A03/P0301240 P040124001401 : No data
/hxmt/work/HXMT-DATA/1L/A03/P0301240 P040124001501 : No data
/hxmt/work/HXMT-DATA/1L/A03/P0301240 P040124001601 : No data
P040124001701
P040124001801
P040124002001
P040124001901
/hxmt/work/HXMT-DATA/1L/A03/P0301240 P040124001701 : No data
P040124002101
/hxmt/work/HXMT-DATA/1L/A03/P0301240 P040124001801 : No data
/hxmt/work/HXMT-DATA/1L/A03/P0301240 P040124002001 : No data
/hxmt/work/HXMT-DATA/1L/A03/P0301240 P040124001901 : No data
/hxmt/work/HXMT-DATA/1L/A03/P0301240 P040124002101 : No data
P040124002201
P040124002301
P040124002401
P040124002501
P040124002601
/hxmt/work/HXMT-DATA/1L/A03/P0301240 P040124002201 : No data
/hxmt/work/HXMT-DATA/1L/A03/P0301240 P040124002301 : No data
/hxmt/work/HXMT-DATA/1L/A03/P0301240 P040124002401 : No data
/hxmt/work/HXMT-DATA/1L/A03/P0301240 P040124002501 : No data
/hxmt/work/HXMT-DATA/1L/A03/P0301240 P040124002601 : No data
P040124002701
P040124002801
P040124002901
P040124003001
P040124003101
/hxmt/work/HXMT-DATA/1L/A03/P0301240 P040124002701 : No data
/hxmt/work/HXMT-DATA/1L/A03/P0301240 P040124002801 : No data
/hxmt/work/HXMT-DATA/1L/A03/P0301240 P040124002901 : No data
/hxmt/work/HXMT-DATA/1L/A03/P0301240 P040124003001 : No data
/hxmt/work/HXMT-DATA/1L/A03/P0301240 P040124003101 : No data
P040124003701
P040124005801
P040124005901
P040124006001
/hxmt/work/HXMT-DATA/1L/A03/P0301240 P040124003701 : No data
P040124006101
/hxmt/work/HXMT-DATA/1L/A03/P0301240 P040124005801 : No data
/hxmt/work/HXMT-DATA/1L/A03/P0301240 P040124005901 : No data
/hxmt/work/HXMT-DATA/1L/A03/P0301240 P040124006001 : No data
/hxmt/work/HXMT-DATA/1L/A03/P0301240 P040124006101 : No data
P040124006201
P040124006301
P040124006401
P040124006501
/hxmt/work/HXMT-DATA/1L/A03/P0301240 P040124006201 : No data
P040124006601
/hxmt/work/HXMT-DATA/1L/A03/P0301240 P040124006301 : No data
/hxmt/work/HXMT-DATA/1L/A03/P0301240 P040124006401 : No data
/hxmt/work/HXMT-DATA/1L/A03/P0301240 P040124006501 : No data
/hxmt/work/HXMT-DATA/1L/A03/P0301240 P040124006601 : No data
P040124006701
P040124006801
P040124006901
P040124007001
/hxmt/work/HXMT-DATA/1L/A03/P0301240 P040124006701 : No data
/hxmt/work/HXMT-DATA/1L/A03/P0301240 P040124006801 : No data
/hxmt/work/HXMT-DATA/1L/A03/P0301240 P040124006901 : No data
/hxmt/work/HXMT-DATA/1L/A03/P0301240 P040124007001 : No data
```