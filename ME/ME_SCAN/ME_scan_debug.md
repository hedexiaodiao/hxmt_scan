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