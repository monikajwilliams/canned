#!/usr/bin/env python
import matplotlib
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import scipy.stats as stats
import numpy as np
import os, sys, math


def smooth_fft(
               positions,
               ):

    fs2ps = 0.001
    stepsize = 50.0
    frate = 1.0/(stepsize*fs2ps)
    sp = np.fft.fft(positions)
    freq = np.fft.fftfreq(positions.shape[-1],frate)

    transfer_ps = 7.0
    transfer_step = transfer_ps*fs2ps*stepsize

    x = 0.00002
#    x = 0.00003
#    x = 0.00004
#    x = 0.00005

    sp2 = np.copy(sp)
    freq2 = np.copy(freq)

    for i, val in enumerate(freq):
        if abs(val) > x:
            sp2[i] = 0.0

    smooth = np.real(np.fft.ifft(sp2))

    return smooth

def test():
    
    fn_in = 'Ostar'
    fs2ps = 0.001
    stepsize=50.0

    ostar = np.load(fn_in)
    tmax = len(ostar)*stepsize*fs2ps
    t = np.linspace(0.0,tmax,len(ostar))
    smoothed_Z = smooth_fft(ostar[:,2])


    dts = np.arange(1,50000,200)
    stds = np.zeros_like(dts)

    plt.clf()
    plt.figure()
    val = 0.003
    for c_count,dt in enumerate(dts):
        run = stepsize*fs2ps*dt
        velocity = np.array([(smoothed_Z[ind+dt]-smoothed_Z[ind])/run for ind in range(len(smoothed_Z)-dt)])
        t = np.linspace(0.0,tmax,len(velocity))
        stds[c_count] += np.std(abs(velocity))
        plt.plot(t,velocity,c=[1.0-c_count*val,0.0,0.0+c_count*val])
        #if c_count != 1:
        #    plt.plot(t,velocity,c=[1.0-c_count*val,0.0,0.0+c_count*val])
        #else:
        #    plt.plot(t,velocity,c='g',linewidth=2)

    plt.axis([0.0,5000.0,-0.3,0.3])
    plt.savefig('velocity_window.pdf')



test()
