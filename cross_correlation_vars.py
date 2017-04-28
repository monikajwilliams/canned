#!/usr/bin/env python
import matplotlib
import math
import numpy as np
import scipy.signal
import matplotlib.pyplot as plt
import pickle

def defects(fn_d = "L_defect",):

    # => Variables <= #    
   
    # Filenames of Data to Load 
    fn_O = "Ostar"
    
   
    # => Loading xyz data for L & D defects <= # 

    print "Loading Data..."
    defect =  pickle.load(open(fn_d,"rb"))
    Ostar = np.load(fn_O)
    print "Data Loaded"

    # => Declaring Lists <= #
   
    # Number of Defects in each Frame 
    num_d = []

    # Z Position of nearest defects to Ostar
    pos_d = []

    # Distance of nearest defects to Ostar
    dist_d = []

    # Z Postions of all defects 
    all_d = []

    # => Looping through Timesteps <= #
    
    for step in range(len(defect)):
        
        # If there are defects in a frame
        if np.shape(defect[step][:])[0] > 0:
            all_Ds = np.array(defect[step][:])[:,2]
            dist_O = abs(all_Ds - Ostar[step,2])

            all_d.append(all_Ds)

            position = all_Ds[np.argsort(dist_O)[0]]
            pos_d.append(position)
            dist_d.append(np.sort(dist_O)[0])
    
            defect_count = len(all_Ds)
            num_d.append(defect_count)
    
        # If there is no defect in a frame
        else:
            num_d.append(0)
            pos_d.append(0)
            dist_d.append(0)
            all_d.append(0)
    
    return all_d, pos_d, dist_d, num_d 

def smooth_fft(
               data,
               # Frequency to Eliminate
               x = 0.00002,
               stepsize = 50.0, # fs 
               visual_verify = False,
               fn_out = 'visual_verify.png',
               ):
    

    # => Smoothing Parameters <= #

    # Units    
    fs2ps = 0.001
    frate = stepsize*fs2ps

    # => Fourier Transform <= #

    # Symmetrizing Data
    #data = np.hstack([data,data[::-1]]) - data.mean()
    data = np.hstack([data,data[::-1]])

    # Real Discrete Fourier Transform of Data
    sp = np.fft.rfft(data)
    # Obtaining Frequencies
    freq = np.fft.rfftfreq(data.shape[0],frate)

    # => Eliminating High Frequency Data Points <= #

    sp2 = np.copy(sp)
    freq2 = np.copy(freq)
    
    for i, val in enumerate(freq):
        if abs(val) > x:
            sp2[i] = 0.0

    # => Reverse Fourier Transform <= #

    smooth_data = np.real(np.fft.irfft(sp2))
    smooth_data_2 = smooth_data[:len(smooth_data)/2+1]

    # => Visually Verifying Data Smoothing <= # 

    if visual_verify == True:
        plot_data = data[:len(data)/2+1]
        end_sim = len(plot_data)*frate
        x_axis1 = np.linspace(0.0,end_sim-1,len(smooth_data_2))
        x_axis2 = np.linspace(0.0,end_sim,len(plot_data))

        plt.clf()
        plt.figure()
        plt.plot(x_axis2,plot_data,'b',alpha=0.3)
        plt.plot(x_axis1,smooth_data_2,'g',alpha=1.0)
        plt.show()


        #plt.figure()
        #plt.plot(freq,sp.real,label='real')
        #plt.plot(freq,sp.imag,label='imag')
        ##plt.plot(freq2,sp2,'r',alpha=0.5,label='cut')
        #plt.legend()
        #plt.show()
    return smooth_data,freq,sp

def flip_rate():

    timestep = 50.0 #fs/step 
    fs2ps = 0.001
    nm2a = 10.0
    
    fn_flips = "flips"
    flips = np.load(fn_flips)

    return flips


def cross_correlation(data_1,
                      data_2,
                      stride,
                      fn_out,
                      timestep = 50.0,     # fs/step
                      len_traj_ps = 100.0,    
                      len_interval_ps = 200.0, 
                      plot_corrs = True,
                      ):
  
    

    # => Unit Conversions <= #

    fs2ps = 0.001 
    len_traj = int(len_traj_ps/(fs2ps*timestep))
    len_interval = int(len_interval_ps/(fs2ps*timestep))

    # => Setting up Indices <= #

    ndata = len(data_1)-len_traj
    npoints = ndata/len_interval
    points = np.linspace(0.0,npoints*len_interval,npoints,dtype=int)

    count = 0.0
    all_corrs = []
    nsamples = (len(data_1)/npoints)*(npoints-1)

    taus = np.arange(0,len_traj,stride,dtype=int)
    ntau = len(taus)
    taus_ps = taus*timestep*fs2ps

    if plot_corrs == True:
        plt.clf()
        plt.figure() 

    avg_v1v2s = []
    for point in points:

        count += 1.0
        
        var_1 = data_1[point:point+len_traj]
        var_2 = data_2[point:point+len_traj]

        avg_v1 = np.average(var_1)
        avg_v2 = np.average(var_2)
        avg_v1v2 = np.sum(np.multiply(var_1,var_2))
        avg_v1v2s.append(avg_v1v2)
        
        correlations = []
        for ind,tau in enumerate(taus):
            if tau == 0.0:
                conj_f = np.conj(var_1)
            else:
                conj_f = np.conj(var_1)[:-tau]
            g = var_2[tau:]
            gf = conj_f*g
            corr = np.sum(gf)
            correlations.append(corr)
            print "step: %d/%d\r" % (ind,ntau),

        correlations = np.array(correlations)
        correlations -= (avg_v1*avg_v2)
        correlations /= (avg_v1v2 - (avg_v1*avg_v2))
    
        color = npoints/10.0    
        if plot_corrs == True:
            plt.plot(taus_ps,correlations,'k',alpha=color/npoints)

        if point == 0.0:
            all_corrs = correlations
        else:
            all_corrs += correlations
        print "\ncount: %d/%d\n" % (count,npoints)

    all_corrs = np.array(all_corrs)
    all_corrs /= count

    y=np.zeros_like(taus_ps)

    if plot_corrs == True:

        plt.plot(taus_ps,all_corrs,'r')
        plt.plot(taus_ps,y,'k--',alpha=0.2)
        plt.xlabel(r'$\mathrm{\tau\ (ps)}$',fontsize=18)
        plt.ylabel(r'$\mathrm{<x(t_{0})\ y(t_{0} + \tau)>}$',fontsize=18)
        plt.axis([0.0,len_traj_ps,-1.0,1.0])
        plt.savefig(fn_out)

    return all_corrs

def autocorrelate(data):
    """Full autocorrelation function of `data`."""
    
    assert len(data.shape) == 1
    N = len(data)
    norm = N - np.abs(np.arange(-N+1, N), dtype=float)

    return scipy.signal.fftconvolve(data, data[::-1], mode="full") / norm

#def test():
#
#    timestep = 50.0 #fs/step 
#    fs2ps = 0.001
#    nm2a = 10.0
#
#    fn_O = "Ostar"
#    dipoles = np.load("Dipoles")
#    fields = smooth_fft(np.load("Fields"))
#    ostar = smooth_fft(np.load(fn_O)[:,2])
#    t = np.arange(0.0,len(ostar))*50.0*fs2ps
#    
#    dt = 1
#    run = timestep*fs2ps*dt
#
#    velocity = np.array([(ostar[ind+dt]-ostar[ind])/run for ind in range(len(ostar)-dt)])
#    #rough_fields = np.load("Fields")
#    #rough_ostar = np.load(fn_O)[:,2]
#    #rough_velocity = np.array([(rough_ostar[ind+dt]-rough_ostar[ind])/run for ind in range(len(rough_ostar)-dt)])
#
#    plt.clf()
#    plt.subplot(211)
#    #plt.plot(t[:-1],rough_velocity,alpha=0.2)
#    plt.plot(t[:-1],velocity)
#    plt.axis([0.0,5000.0,-0.5,0.5])
#    plt.ylabel('velocity')
#
#    plt.subplot(212)
#    #plt.plot(t,rough_fields,alpha = 0.2)
#    plt.plot(t,fields)
#    plt.ylabel('E field')
#    plt.xlabel('t (ps)')
#    plt.axis([0.0,5000.0,-10.0,5.0])
#    plt.savefig("vt_vf.png")
#
#    stride = 1
#
#    cross_correlation(fields[:-1],
#                      velocity,
#                      fn_out='fv.png',
#                      len_traj_ps=50,
#                      len_interval_ps=60,
#                      stride=stride
#                      )
##    cross_correlation(fields,
##                      fields,
##                      fn_out='ff.pdf',
##                      stride=stride,
##                      len_traj_ps=100,
##                      len_interval_ps=110,
##                      )  
##    cross_correlation(velocity,
##                      velocity,
##                      fn_out='vv.pdf',
##                      len_traj_ps=200,
##                      len_interval_ps=250,
##                      stride=stride)
#
##    cross_correlation(
##                      velocity,
##                      fields[:-1],
##                      fn_out='vf.pdf',
##                      len_traj_ps=50,
##                      len_interval_ps=60,
##                      stride=stride
##                      )
#
##    random_1 = np.random.rand(10000)
##    random_2 = np.random.rand(10000)
##    line_1 = np.linspace(0.0,10000.0,10000)
##    line_2 = np.linspace(0.0,1.0,10000)
##    test_1 = np.sin(line_1)
##    test_2 = line_1*-1.0
##
##    timestep = 50.0
##    fs2ps = 0.001
##    len_traj_ps = 75
##    len_interval_ps = 100
##    len_traj = len_traj_ps/(fs2ps*timestep)
##    len_interval = len_interval_ps/(fs2ps*timestep)
##
##    plt.clf()
##    plt.figure()
##    plt.plot(line_1[:len_traj]*fs2ps*timestep,test_1[:len_traj],alpha=0.3)
##    plt.plot(line_1[:len_traj]*fs2ps*timestep,test_2[:len_traj],alpha=0.3)
##    plt.savefig('test_plot.png')
##
##    cross_correlation(
##                      test_1,
##                      test_2,
##                      fn_out='test.png',
##                      len_traj_ps=len_traj_ps,
##                      len_interval_ps=len_interval_ps,
##                      stride=stride
##                      )
# 
#
#test()
