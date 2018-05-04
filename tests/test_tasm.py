# -*- coding: utf-8 -*-
"""
Created on Mon Mar 26 15:55:45 2018

@author: Dimitris Loukrezis

Test DALI algorithm on test function
f(x) = exp[-x1^2]*cos(x2)*exp[-x3^2]*cos[x4]
"""

import chaospy as cp
import numpy as np
import sys
sys.path.append("../dali")
from dali import dali
from postprocess import compute_cv_error, compute_moments
import math

def qoi_tasm(x):
    """Meromorphic function"""
    return math.exp(- x[0] ** 2) * math.cos(x[1]) 
#* math.exp(- x[2] ** 2) * math.cos(x[3])


# num RVs
N = 2
# joint PDF
jpdf = cp.Iid(cp.Uniform(-1,1), N)
# function to be approximated
f = qoi_tasm
# maximum function calls
max_fcalls = np.linspace(10, 90, 9).tolist()
max_fcalls = max_fcalls + np.linspace(100, 200, 10).tolist()

# arrays for results storage
cv_errz = []
meanz = []
varz = []
fcallz = []
ratioz = []
# approximate
interp_dict = {}
for mfc in max_fcalls:
    mfc = int(mfc)
    interp_dict = dali(f, N, jpdf, tol=1e-16, max_fcalls=mfc, 
                       interp_dict=interp_dict, verbose=False)
    # cross-validation error
    cv_err = compute_cv_error(interp_dict, jpdf, f, Nsamples=1000)
    cv_errz.append(cv_err)
    mu, sigma2 = compute_moments(interp_dict, jpdf)
    meanz.append(mu)
    varz.append(sigma2)
    fcalls = len(interp_dict['idx_act'] + interp_dict['idx_adm'])
    fcallz.append(fcalls)
    ratio = float(len(interp_dict['idx_act'])) / float(fcalls)
    ratioz.append(ratio)
    print ""
    print "Max. fcalls = ", mfc
    print "Cross-validation error = ", cv_err
    print "Expected value = ", mu
    print "Variance = ", sigma2
    print "Function calls = ", fcalls
    print "Act/Adm ratio = ", ratio
    print ""

print "This is the end, beautiful friend."
