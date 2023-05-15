"""
Search for  a burst in the data of two detectors and calculate cross-correlation lag.
"""

import os
import sys
import numpy as np

from scipy.stats import chi2, norm

def qw(s):
    return s.split()

class cc_parameters():
    """A simple class to store information for cross correlation."""
    def __init__(self):
        self.info = {} 
        self.keys = qw('sc1 sc2 th1 th2 i_beg_2 t_beg_2 i_end_2 t_end_2')

def get_dTcc(arr_dt, arr_chi, nDOF, fRijMin, nMin, nSigma=3):
    """Calculate cross-correlation lag confidence interval

    Parameters
    ----------
    arr_dt 
    arr_chi 
    nDOF 
    fRijMin 
    nMin 
    nSigma
    
    Returns
    -------
    dTlower, dTupper - lower and upper boundary of the interval
    fSigma - chi2/dof level for nSigma
    """
    
    #fSigma = mychi2(nDOF-1, nSigma) - 1.0 + fRijMin

    p_0 = norm.cdf(nSigma) - norm.cdf(-nSigma)
    fSigma = chi2.ppf(p_0, nDOF-1)/(nDOF -1) - 1.0 + fRijMin

    fSigmaSimple = fRijMin + nSigma**2/nDOF # d chi2 = 9 for 3 sigma

    n = arr_chi.size

    #search upper dTcc 3 sigma
    i = nMin
    fRij = arr_chi[i]

    while((i > 1) and (fRij < fSigma)):
        i = i-1
        fRij = arr_chi[i]
    
    dTupper = arr_dt[i]
    
    #search lower dTcc 3 sigma
    i = nMin
    fRij = arr_chi[i]
    while((i < n-1) and (fRij < fSigma)):
        i = i + 1 
        fRij = arr_chi[i]
    
    dTlower = arr_dt[i]
    
    return dTlower, dTupper, fSigma


def correlate(
    arr_t_1, arr_cnts_1, arr_cnts_err_1, 
    arr_t_2, arr_cnts_2, arr_cnts_err_2,
    i_beg_2, i_end_2,        
    i_beg_1, n_max_1,        
    fscale,         
    dt_1,           
    dt_2,           
    n_sum_2=1       
):
    """ Calculate chi2/dof(dTcc)

    Parameters
    ----------
    arr_t_1         array with TH 1 bin starts
    arr_cnts_1      array with TH 1 bg sub counts
    arr_cnts_err_1  array with TH 1 bg sub counts err 
    arr_t_2         array with TH 2 bin starts
    arr_cnts_2      array with TH 2 bg sub counts
    arr_cnts_err_2  array with TH 2 bg sub counts err 
    i_beg_2         start bin of TH 2 to cross-correlate
    i_end_2         end bin of TH 2 to cross-correlate
    i_beg_1         start bin of TH 1
    n_max_1         maximum shift of start of TH 1
    fscale          scale TH 2 to TH 1
    dt_1            time resolution of sc1, ms!
    dt_2            time resolution of sc2, ms!
    n_sum_2=1       number of bins of TH 2 to sum, default = 1

    Returns
    -------
    arr_dt, arr_chi, nDOF, fRijMin, dTmin, iMin, nMin

    Notes
    -----
    We shift the time history 1 !!!!
    the time history 2 is in rest
    """

    arr_dt = np.zeros(n_max_1)
    arr_chi = np.zeros(n_max_1)
    arr_sigma3 = np.zeros(n_max_1)
    arr_sigma2 = np.zeros(n_max_1)
    arr_sigma3simple = np.zeros(n_max_1)
    arr_sigma2simple = np.zeros(n_max_1)
    
    n = 0   
    fRijMin = 1e4
    i_min = 0
    n_min = 0

    nDOF = (i_end_2 - i_beg_2) / n_sum_2 # number degrees of freedom
    fRij = 0

    for i in range(i_beg_1, i_beg_1 + n_max_1):
    
        fRij = 0
    
        fT1 = dt_1
        fT2 = 0
        
        l = i
        
        fdCnt1 = 0
        fdErr1 = 0
        
        for k in range(i_beg_2, i_end_2 + 1, n_sum_2):
        
            fCnt2 = 0
            fErr2 = 0
            
            for iSum2 in range(k, k + n_sum_2):
                fCnt2 += arr_cnts_2[iSum2]
                fErr2 += arr_cnts_err_2[iSum2]
            
            fCnt2 *= fscale
            fErr2 *= fscale * fscale
            
            fCnt1 = fdCnt1
            fErr1 = fdErr1
            
            while (fT1 <= (fT2 + dt_2 * n_sum_2)):
                fCnt1 += arr_cnts_1[l]
                fErr1 += arr_cnts_err_1[l]
                fT1 += dt_1
                l = l + 1
            
            fdT = (fT1 - (fT2 + dt_2 * n_sum_2)) / dt_1
            fdCnt1 = fdT * arr_cnts_1[l]
            fCnt1 += (1.0 - fdT) * arr_cnts_1[l]
            
            fdErr1 = fdT * arr_cnts_err_1[l]
            fErr1 += (1.0 - fdT) * arr_cnts_err_1[l]
            
            l = l + 1
            fT1 += dt_1
            fT2 += dt_2 * n_sum_2         
            
            fDif = fCnt1 - fCnt2
            fDif *= fDif
            
            fDif /= (fErr1 + fErr2)
            fRij += fDif
       
        fRij /= (nDOF - 1)
        
        dT = arr_t_2[i_beg_2] - arr_t_1[i]
        arr_dt[n] = dT
        arr_chi[n] = fRij
        
        if (fRij < fRijMin):
            fRijMin = fRij
            iMin = i
            nMin = n
            dTmin = dT  
        
        n = n + 1
    
    return arr_dt, arr_chi, nDOF, fRijMin, dTmin, iMin, nMin

def out_x_y(x,y, file_name):
   """
   A simple test output of two related numpy arrays
   """
   with open(file_name, 'w') as f:
       for i in range(x.size):
          f.write("{:8.3f} {:8.3f}\n".format(x[i], y[i]))
   
def test_array(a, desc=None):
    """
    A simple test output of a numpy array with a description
    """
    np.set_printoptions(suppress=True)
    if desc:
        print (desc, a, a.size)
    else:
        print (a, a.size)

 