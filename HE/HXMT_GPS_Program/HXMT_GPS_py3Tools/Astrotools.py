#!/hxmt/soft/Develop/anaconda2/bin/python
import numpy as np
import math
#import matplotlib
#matplotlib.use('Agg')
#matplotlib.rcParams['xtick.direction'] = 'in'
#matplotlib.rcParams['ytick.direction'] = 'in'
#import matplotlib.pyplot as plt
#from astropy.io import fits as pf
#import sys,os
#from glob import glob as glob

def distance(ratpp,dectpp,rapnt,decpnt):
    ytp = (np.sin(decpnt*np.pi/360.0 - dectpp*np.pi/360.0))**2 + np.cos(decpnt*np.pi/180.0)*np.cos(dectpp*np.pi/180.0)*(np.sin(ratpp*np.pi/360.0 - rapnt*np.pi/360.0))**2
    ytp = np.arccos(1-2*ytp)*180./np.pi
    return ytp


class Measurement(object):
    def __init__(self, val, perc):
        self.val = val
        self.perc = perc
        self.abs  = self.val * self.perc / 100.0
    def __repr__(self):
        return "Measurement (%r, %r)" % (self.val, self.perc)
    def __str__(self):
        return "%g+-%g%%" % (self.val, self.perc)
    def _addition_result(self, result, other_abs):
        new_perc = 100.0 *(math.hypot(self.abs, other_abs) / result)
        return Measurement(result, new_perc)
    def __add__(self, other):
        result = self.val + other.val
        return self._addition_result(result, other.abs)
    def __sub__(self, other):
        result = self.val - other.val
        return self._addition_result( result, other.abs)
    def _multiplication_result(self, result, other_perc):
        new_perc = math.hypot(self.perc, other_perc)
        return Measurement(result, new_perc)
    def __mul__(self, other):
        result = self.val * other.val
        return self._multiplication_result(result, other.perc)
    def __div__(self, other):
        result = self.val / other.val
        return self._multiplication_result(result, other.perc)

def Sel_Median(array,sigma = 3.):
    x = np.copy(array)
    od = np.arange(len(x))
    while 1:
        md = np.median(x)
        td = np.std(x)
        bad = x > md + sigma*td
        if bad.sum()>0:
            x = x[~bad]
            od = od[~bad]
        else:
            return x,od


