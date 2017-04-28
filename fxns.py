#!/usr/bin/env python
import matplotlib
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import scipy.stats as stats
import numpy as np
import os, sys, math


def cauchy(x,x_0,gamma,height):

    a = 1.0/(math.pi * gamma)
    b = (x-x_0)**2
    c = gamma**2
    norm = max(a*(c/(b+c)))

    return a*(c/(b+c))*(height/norm)

def gaussian(x,mu,sigma,height):

    a = 1.0/np.sqrt(2.0*math.pi*(sigma**2))
    b = np.exp(-((x-mu)**2)/(2.0*(sigma**2)))
    norm = max(a*b)

    return a*b*(height/norm)

def exp_1(x,A,T):

    val = A*np.exp(-x/T)

    return val
    
def exp_2(x,A_1,T_1,A_2,T_2):

    val = A_1*np.exp(-x/T_1) + A_2*np.exp(-x/T_2)

    return val

def const_exp_(x,A_1,T_1,A_2):

    val = A_1*np.exp(-x/T_1) + A_2

    return val
    
def exp_2_osc(x,A_1,T_1,A_2,T_2,A_3,T_3,omega):

    val1 = A_1*np.exp(-x/T_1) + A_2*np.exp(-x/T_2) 
    val2 = A_3*np.cos(omega*x)
    val3 = (1.0/(T_3*omega))*np.sin(omega*x)*np.exp(-x/T_3)
    val = val1 + val2 + val3

    return val

def exp_2_osc_const(x,A_1,A_2,T_2,A_3,T_3,omega):

    val1 = A_1 + A_2*np.exp(-x/T_2) 
    val2 = A_3*np.cos(omega*x)
    val3 = (1.0/(T_3*omega))*np.sin(omega*x)*np.exp(-x/T_3)
    val = val1 + val2 + val3

    return val

def fit_func(func,
             xdata,
             ydata,
             p0,
             ):
    
    popt,pcov = curve_fit(func,xdata,ydata,p0=p0)
    perr = np.sqrt(np.diag(pcov))

    #print "Standard Deviation Error on Parameters:"
    #for p in perr:
    #    print "%5.5E" % (p)

    return popt,perr
    
