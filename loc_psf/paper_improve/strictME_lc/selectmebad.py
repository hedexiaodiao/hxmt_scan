from astropy.io import fits as pf
import numpy as np
from scipy.optimize import curve_fit

def gaussian(x,*param):
    return param[0]*np.exp(-np.power(x - param[1], 2.) / (2 * np.power(param[2], 2.)))


def gaussfit(x,y):
    popt,pcov = curve_fit(gaussian,x,y)
    sigma = np.sqrt(np.diag(pcov))
    y1=gaussian(x,*popt)
    return y1

