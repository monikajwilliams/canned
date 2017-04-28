#!/usr/bin/env python
import os, sys, re, math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from canned import cross_correlation_vars as ccv
from canned import fxns
from traj_analysis import trackOstar as trkO
import time
from datetime import date

def plot_correlations(
                     fn_data_1,
                     fn_data_2,
                     fn_out,
                     smooth_1 = False,
                     smooth_2 = False,
                     velocity_1 = False,
                     velocity_2 = False,
                     select = 10000000,
                     timestep = 50.0,
                     len_traj_ps = 200,
                     len_interval = 100,
                     dt = 1,
                     stride = 1,
                     plot_axis = [0.0],
                     axis = 2,
                     ):
    densities = [
            #3.0,
            #3.2,
            #3.4,
            #3.5,    
            #3.6,
            #3.7,
            #3.8,
            3.9,
            #4.0,
            #4.1,
            #4.2,
            #4.3,
            #4.4,
            #4.5,
             ]
        
    
    samples = [
               0,
               1,
               2,
               3,
               4,
               5,
               6,
               7,
               8,      
               9
                ]

    date_1 = date.today()
    timestamp = "%s-%s" % (str(date_1.day), str(date_1.month))
    ncount = 0
    
    nd = 1
    fs2ps = 0.001
    nm2a = 10.0
    run = timestep*fs2ps*dt
    
    len_traj = len_traj_ps/(fs2ps*timestep)
    #len_interval = len_interval_ps/(fs2ps*timestep)
    len_interval_ps = len_interval*(fs2ps*timestep)
    taus = np.arange(0,len_traj,stride)
    ntau = len(taus)
    taus_ps = taus*timestep*fs2ps
    
    all_corrs = []
    for m in densities:
        for n in samples:
    
            if m == 4.5 and n+1 <= 9:
                continue
            elif m == 4.4 and n+1 == 1:
                continue
            elif m == 3.0 and n+1 <= 9:
                continue
            elif m == 3.5 and n+1 == 1:
                continue
            elif m == 3.7 and n+1 == 1:
                continue
            elif m == 3.7 and n+1 == 8:
                continue
            elif m == 3.8 and n+1 == 5:
                continue
            else:
    
    	        os.chdir('/Users/monikawilliams/Desktop/current_Research/stats_100/RAN_stats_100_%1.1fD/PD%d/' % (m,n+1))
    	        os.system("pwd")
                fn_O = 'Ostar'

                print "Loading Data"
                if fn_data_1 == fn_O:
                    data_1 = np.load(fn_data_1)[:,axis]
                    ostar = data_1
                    t = np.arange(0.0,len(ostar))*timestep*fs2ps
                else:
                    data_1 = np.load(fn_data_1)

                print "Loading Data"
                if fn_data_2 == fn_O:
                    data_2 = np.load(fn_data_2)[:,axis]
                    t = np.arange(0.0,len(ostar))*timestep*fs2ps
                else:
                    data_2 = np.load(fn_data_2) 

                if smooth_1 == True:
                    print "Smoothing Data"
                    if fn_data_1 == fn_O:
                        ostar = ccv.smooth_fft(ostar)
                    else:
                        data_1 = ccv.smooth_fft(data_1)
                if smooth_2 == True:
                    print "Smoothing Data"
                    if fn_data_2 == fn_O and fn_data_1 != fn_O:
                        ostar = ccv.smooth_fft(ostar)
                    else:
                        data_2 = ccv.smooth_fft(data_2)

                
                if velocity_1 == True:                
                    data_1 = np.array([(ostar[ind+dt]-ostar[ind])/run for ind in range(len(ostar)-dt)])
                    if velocity_2 != True:
                        data_2 = data_2[:-1]
                if velocity_2 == True:                
                    data_2 = np.array([(ostar[ind+dt]-ostar[ind])/run for ind in range(len(ostar)-dt)])
                    if velocity_1 != True:
                        data_1 = data_1[:-1]
    

                sample_cross = ccv.cross_correlation(
                                  data_1,
                                  data_2,
                                  len_traj_ps=len_traj_ps,
                                  len_interval_ps=len_interval_ps,
                                  fn_out=fn_out,
                                  stride=stride,
                                  plot_corrs = False,
                                  )
                if np.max(abs(sample_cross)) > select:
                    print "GOTCHA"
                    os.system("pwd")
                else:
                    if ncount == 0:
                        all_corrs = sample_cross
                    else:
                        all_corrs += sample_cross
                    plt.plot(taus_ps,sample_cross,'k',alpha=0.1)
                    ncount += 1
    
    os.chdir('/Users/monikawilliams/Desktop/current_Research/stats_100/scripts/plots/')
    y=np.zeros_like(taus_ps)
    all_corrs /= ncount
    plt.plot(taus_ps,all_corrs,'r')
    plt.plot(taus_ps,y,'k--',alpha=0.2)
    plt.xlabel(r'$\mathrm{\tau\ (ps)}$',fontsize=18)
    #plt.ylabel(r'$\mathrm{<x(t_{0})\ y(t_{0} + \tau)>}$',fontsize=18)
    plt.ylabel(r'$\mathrm{C(t)}$',fontsize=18)
    if len(plot_axis) > 1:
        plt.axis(plot_axis)
    plt.savefig(fn_out)


    return all_corrs


def test():

    len_traj_ps=50
    len_interval=1000
    timestep = 50.0
    fs2ps = 0.001
    stride = 1
    len_traj = len_traj_ps/(fs2ps*timestep)
    len_interval_ps = len_interval*(fs2ps*timestep)
    taus = np.arange(0,len_traj,stride)
    ntau = len(taus)
    taus_ps = taus*timestep*fs2ps

    vv_corrs =  plot_correlations(
                      fn_data_1='fields',
                      fn_data_2='fields',
                      fn_out='fit_rr_ff_1000.png',
                      smooth_1 = False,
                      smooth_2 = False,
                      velocity_1 = False,
                      velocity_2 = False,
                      select = 10000000,
                      timestep = 50.0,
                      len_traj_ps = len_traj_ps,
                      len_interval = len_interval,
                      dt = 1,
                      stride = 1,
                      plot_axis = [0.0,len_traj_ps,-1.0,1.0],
                      axis = 2,
                      )

    guess = [1.0,50.0,1.0,20.0]
    vals,err = fxns.fit_func(fxns.exp_2,taus_ps,vv_corrs,p0=guess)
    fit_vals = fxns.exp_2(taus_ps,vals[0],vals[1],vals[2],vals[3]) 
    print "fit 2 exp"
    print vals

    guess_1 = [1.0,50.0]
    vals_1,err_1 = fxns.fit_func(fxns.exp_1,taus_ps,vv_corrs,p0=guess_1)
    fit_vals_1 = fxns.exp_1(taus_ps,vals_1[0],vals_1[1]) 
    print "fit 1 exp"
    print vals_1

    guess_2 = [1.0,1.0,1.0,0.2,1.0,0.2,30.0]
    vals_2,err_2 = fxns.fit_func(fxns.exp_2_osc,taus_ps,vv_corrs,p0=guess_2)
    fit_vals_2 = fxns.exp_2_osc(taus_ps,vals_2[0],vals_2[1],vals_2[2],vals_2[3],
                                vals_2[4],vals_2[5],vals_2[6]) 
    print "fit 2 exp osc"
    print vals_2

    guess_3 = [1.0,1.0,0.2,1.0,0.2,30.0]
    vals_3,err_3 = fxns.fit_func(fxns.exp_2_osc_const,taus_ps,vv_corrs,p0=guess_3)
    fit_vals_3 = fxns.exp_2_osc_const(taus_ps,vals_3[0],vals_3[1],vals_3[2],vals_3[3],
                                vals_3[4],vals_3[5]) 
    print "fit 2 exp osc const"
    print vals_3

    plt.clf()
    plt.figure()
    plt.plot(taus_ps,vv_corrs,label=r'$\mathrm{E_{z}}$')
    plt.plot(taus_ps,fit_vals,label=r'$\mathrm{A_{1}e^{\frac{-t}{\tau_{1}}} + A_{2}e^{\frac{-t}{\tau_{2}}}}$')
    plt.plot(taus_ps,fit_vals_1,label=r'$\mathrm{Ae^{\frac{-t}{\tau}}}$')
    #plt.plot(taus_ps,fit_vals_2,label=r'$\mathrm{fit 2 exp osc}$')
    #plt.plot(taus_ps,fit_vals_3,label=r'$\mathrm{fit 2 exp osc const}$')
    #plt.legend(fontsize=16)
    plt.axis([0.0,2.0,-1.0,1.0])
    plt.savefig('zoom_smooth_fit_ff_%f.ff.pdf' % (len_traj_ps))


   #vv_corrs = plot_correlations(
   #                   fn_data_1='Ostar',
   #                   fn_data_2='Ostar',
   #                   fn_out='rough_vv_10.png',
   #                   smooth_1 = False,
   #                   smooth_2 = False,
   #                   velocity_1 = True,
   #                   velocity_2 = True,
   #                   select = 10000000,
   #                   timestep = 50.0,
   #                   len_traj_ps = len_traj_ps,
   #                   len_interval = len_interval,
   #                   dt = 1,
   #                   stride = 1,
   #                   plot_axis = [0.0,len_traj_ps,-1.0,1.0],
   #                   axis = 2,
   #                   )

    
        
    

test()
