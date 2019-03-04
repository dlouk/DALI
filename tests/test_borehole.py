# -*- coding: utf-8 -*-
"""
Created on Thu Mar 22 11:18:16 2018

@author: Dimitris Loukrezis

Test the DALI algorithm on a borehole model.
"""

import chaospy as cp
import numpy as np
import sys
sys.path.append("../dali")
from dali import dali
from postprocess import compute_cv_error, compute_moments


def borehole(x):
    """Analytical model of the water flow through a borehole"""
    rw = x[0] # radius of borehole (m)
    r  = x[1] # radius of influence (m)
    Tu = x[2] # transmissivity of upper aquifer (m^2/yr)
    Hu = x[3] # potentiometric head of upper aquifer (m)
    Tl = x[4] # transmissivity of lower aquifer (m^2/yr)
    Hl = x[5] # potentiometric head of lower aquifer (m)
    L  = x[6] # length of borehole (m)
    Kw = x[7] # hydraulic conductivity of borehole (m/yr)
    #
    frac_nom = 2 * np.pi * Tu * (Hu - Hl)
    frac_help = 2*L*Tu/(np.log(r/rw)*rw*rw*Kw)
    frac_denom = np.log(r/rw) * (1 + frac_help + Tu/Tl)
    #
    response = frac_nom / frac_denom
    return response
#bot_inp = [0.05, 100., 63070., 990., 63.1, 700., 1120., 9855.]
#print borehole(bot_inp)
#top_inp = [0.15, 50000., 115600., 1110., 116., 820., 1680., 12045.]
#print borehole(top_inp)

# parameter PDFs
# cp.Truncnorm(lo, up, mu,sigma)
pdf1 = cp.Truncnorm(0.05, 0.15, 0.1, 0.0165)
pdf2 = cp.Truncnorm(100, 50000, 3700, 4890)
pdf3 = cp.Truncnorm(63070., 115600.)
pdf4 = cp.Truncnorm(990., 1110., 1050., 20.)
pdf5 = cp.Truncnorm(63.1, 116., 90., 9.)
pdf6 = cp.Truncnorm(700., 820., 760., 20.)
pdf7 = cp.Truncnorm(1120., 1680., 1400., 100.)
pdf8 = cp.Truncnorm(9855., 12045., 10950., 400.)
#
# joint PDF
jpdf = cp.J(pdf1, pdf2, pdf3, pdf4, pdf5, pdf6, pdf7, pdf8)
# num RVs
N = 8
## function to be approximated
f = borehole
# results storage
max_fcalls = np.linspace(10,90,9).tolist()
max_fcalls = max_fcalls + np.linspace(100, 1000, 10).tolist()
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
                       interp_dict=interp_dict, verbose=True)
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