#!/usr/bin/env python
import os, sys, re, math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from canned import cross_correlation_vars as ccv
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
                     acceleration_1 = False,
                     acceleration_2 = False,
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
               #1,
               #2,
               #3,
               #4,
               #5,
               #6,
               #7,
               #8,      
               #9
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
    
    	        #os.chdir('/Users/monikawilliams/Desktop/current_Research/stats_100/RAN_stats_100_%1.1fD/PD%d/' % (m,n+1))
    	        os.system("pwd")
                fn_O = 'rs'

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
                    ostar = data_2
                    t = np.arange(0.0,len(ostar))*timestep*fs2ps
                else:
                    data_2 = np.load(fn_data_2) 


                n_avg = 1
                if smooth_1 == True:
                    print "Smoothing Data"
                    if fn_data_1 == fn_O:
                  
                        print "locally averaging positions" 
                        smooth_positions = np.zeros(np.shape(ostar))     
                        for i,val in enumerate(ostar):
                            if i > n_avg and i < len(ostar)-2.0*n_avg:
                                smooth_positions[i] = np.mean(ostar[i-n_avg:i+n_avg])

                        ostar2 = ccv.smooth_fft(smooth_positions,
                        #ostar2 = ccv.smooth_fft(ostar,
                                               visual_verify=True,
                                               fn_out ='test_fft.png',
                                               x = 0.008,
                                                )
                    else:
                        data_1 = ccv.smooth_fft(data_1)
                if smooth_2 == True:
                    print "Smoothing Data"
                    if fn_data_2 == fn_O and fn_data_1 != fn_O:
                        data_2 = np.array([(ostar[ind+dt]-ostar[ind])/run for ind in range(len(ostar)-dt)])
                        ostar2 = ostar
                        ostar3,freq1,sp1 = ccv.smooth_fft(ostar2,
                        #vel_test,freq1,sp1 = ccv.smooth_fft(data_2,
                                               visual_verify=True,
                                               fn_out ='test_fft.png',
                                               x = 0.007,
                                               )

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

                if acceleration_1 == True:
                    vel_1 = np.array([(ostar[ind+dt]-ostar[ind])/run for ind in range(len(ostar)-dt)])
                    data_1 = np.array([(vel_1[ind+dt]-vel_1[ind])/run for ind in range(len(vel_1)-dt)])
                    if acceleration_2 != True and velocity_2 != True:
                        data_2 = data_2[:-2]
                    elif velocity_2 == True:
                        data_2 = data_2[:-1]

                if acceleration_2 == True:
                    vel_2 = np.array([(ostar[ind+dt]-ostar[ind])/run for ind in range(len(ostar)-dt)])
                    data_2 = np.array([(vel_2[ind+dt]-vel_2[ind])/run for ind in range(len(vel_2)-dt)])
                    if acceleration_1 != True and velocity_1 != True:
                        data_1 = data_1[:-2]
                    elif velocity_1 == True:
                        data_1 = data_1[:-1]
    
                #sample_cross = ccv.cross_correlation(
                #                  data_1,
                #                  data_2,
                #                  len_traj_ps=len_traj_ps,
                #                  len_interval_ps=len_interval_ps,
                #                  fn_out=fn_out,
                #                  stride=stride,
                #                  plot_corrs = False,
                #                  )
                sample_cross = ccv.autocorrelate(data_1)
                t = np.arange(len(sample_cross))
                plt.figure()
                plt.plot(t,sample_cross)

                plt.figure()
                plt.plot(np.fft.rfft(sample_cross))
                plt.show()

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
    
    #os.chdir('/Users/monikawilliams/Desktop/current_Research/stats_100/scripts/plots/')
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


def test():

#    plot_correlations(
#                      fn_data_1='fields',
#                      fn_data_2='fields',
#                      fn_out='rough_ff_1.png',
#                      smooth_1 = False,
#                      smooth_2 = False,
#                      velocity_1 = False,
#                      velocity_2 = False,
#                      select = 10000000,
#                      timestep = 50.0,
#                      len_traj_ps = 1,
#                      len_interval = 1,
#                      dt = 1,
#                      stride = 1,
#                      plot_axis = [0.0,0.6,-1.0,1.0],
#                      axis = 2,
#                      )

#    plot_correlations(
#                      fn_data_1='Ostar',
#                      fn_data_2='Ostar',
#                      fn_out='smooth_vv_10.png',
#                      smooth_1 = True,
#                      smooth_2 = True,
#                      velocity_1 = True,
#                      velocity_2 = True,
#                      select = 10000000,
#                      timestep = 50.0,
#                      len_traj_ps = 10.0,
#                      len_interval = 1,
#                      dt = 1,
#                      stride = 1,
#                      plot_axis = [0.0,10.0,-1.0,1.0],
#                      axis = 2,
#                      )
#
#    plot_correlations(
#                      fn_data_1='fields',
#                      fn_data_2='Ostar',
#                      fn_out='rough_smooth_fv_1.png',
#                      smooth_1 = False,
#                      smooth_2 = True,
#                      velocity_1 = False,
#                      velocity_2 = True,
#                      select = 10000000,
#                      timestep = 50.0,
#                      len_traj_ps = 1.0,
#                      len_interval = 1,
#                      dt = 1,
#                      stride = 1,
#                      plot_axis = [0.0,1.0,-3.0,3.0],
#                      axis = 2,
#                      )
#    plot_correlations(
#                      fn_data_1='Ostar',
#                      fn_data_2='Ostar',
#                      fn_out='rough_aa_5.png',
#                      smooth_1 = False,
#                      smooth_2 = False,
#                      velocity_1 = False,
#                      velocity_2 = False,
#                      acceleration_1 = True,
#                      acceleration_2 = True,
#                      select = 10000000,
#                      timestep = 50.0,
#                      len_traj_ps = 5.0,
#                      len_interval = 100,
#                      dt = 1,
#                      stride = 1,
#                      plot_axis = [0.0],
#                      axis = 2,
#                      )

    plot_correlations(
                      fn_data_1='rs',
                      fn_data_2='rs',
                      fn_out='rr_pp_500.png',
                      smooth_1 = False,
                      smooth_2 = False,
                      velocity_1 = False,
                      velocity_2 = False,
                      acceleration_1 = False,
                      acceleration_2 = False,
                      select = 100000,
                      timestep = 50.0,
                      len_traj_ps = 500.0,
                      len_interval = 1000,
                      dt = 1,
                      stride = 1,
                      plot_axis = [0.0],
                      axis = 2,
                      )
    #plot_correlations(
    #                  fn_data_1='Ostar',
    #                  fn_data_2='Ostar',
    #                  fn_out='ss_aa_1000.png',
    #                  smooth_1 = True,
    #                  smooth_2 = True,
    #                  velocity_1 = False,
    #                  velocity_2 = False,
    #                  acceleration_1 = True,
    #                  acceleration_2 = True,
    #                  select = 10000000,
    #                  timestep = 50.0,
    #                  len_traj_ps = 1000,
    #                  len_interval = 3000,
    #                  dt = 1,
    #                  stride = 1,
    #                  plot_axis = [0.0],
    #                  axis = 2,
    #                  )
test()
