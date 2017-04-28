#!/usr/bin/env python
import matplotlib
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import scipy.stats as stats
import numpy as np
import os, sys, math

def velocities(traj_len = 5000.0,
               ntraj = 11,
               nd = 1,
               axis = 2,
               fn_positions = "Ostar",
               n_avg = 500,
               ):

    nm2a = 10.0
    fn_velocities = "velocities.png"
    simulation_len = traj_len*ntraj #ps
    skip_traj = 10000

    print "Loading position data"
    positions = []
    nsteps = 0
    ntraj = 1
    for n in range(ntraj):
        fn_in = "%s%d" % (fn_positions,n)
        if nd == 3:
            position = np.array(np.load(fn_in))[:,:]*nm2a
            position[:,0] -= min(position[:,0])
            position[:,1] -= min(position[:,1])
            position[:,2] -= min(position[:,2])
        elif nd == 1:
            position = np.array(np.load(fn_in))[:,2]*nm2a
            position -= min(position)
        else:
            "Check Your Dimensionality! - 2D system still in progress"
        step2ps = simulation_len/(len(position)) #ps/steps
        positions.append(position)
        print "Loaded positions %d/%d" % (n,ntraj)
        nsteps = len(position)

    positions = np.reshape(np.array(positions),(nsteps*ntraj,nd))
    print "Smoothing position data"
    if nd == 1:
        smooth_positions = np.zeros(len(positions)-2.0*n_avg)
        times_vel = np.linspace(0.0,simulation_len,len(smooth_positions))
        times_pos = np.linspace(0.0,simulation_len,len(positions))
        for i,val in enumerate(positions):
            if i > n_avg and i < len(positions)-2.0*n_avg:
                smooth_positions[i] = np.mean(positions[i-n_avg:i+n_avg])

        dy = step2ps
        print "Taking gradient of position data"
        velocities = np.gradient(smooth_positions[skip_traj:],dy)
        
        plt.clf()
        plt.figure()
        plt.subplot(212)
        plt.plot(times_vel[skip_traj:],velocities,label="velocity")
        plt.ylabel(r"$\mathrm{Defect\ Velocity\ \AA/ps}$",fontsize=16) 
        plt.xlabel(r"$\mathrm{Time\ ps}$",fontsize=16) 
        #plt.xlabel(r"$\mathrm{Time\ ps}$") 

        plt.subplot(211)
        plt.plot(times_pos,positions,label="positions",alpha=0.2)
        plt.plot(times_vel,smooth_positions,label="smooth positions")
        plt.ylabel(r"$\mathrm{Defect\ Position\ \AA}$",fontsize=16) 
        #plt.legend()
        plt.savefig(fn_velocities)

        plt.clf()
        plt.figure()
        plt.plot(times_pos[:20000], positions[skip_traj:30000],alpha=0.2)
        plt.plot(times_vel[:20000], smooth_positions[skip_traj:30000])
        plt.plot(times_vel[:20000], velocities[:20000])
        plt.savefig("segment_fit.png")

        plt.clf()
        plt.figure()
        indsVel = np.argsort(velocities)
        pos_vel = positions[skip_traj:]
        plt.plot(pos_vel[indsVel], velocities[indsVel],'o',alpha=0.2)
        plt.xlabel(r"$\mathrm{Positions\ \AA}$",fontsize=18)
        plt.ylabel(r"$\mathrm{Velocities\ \AA/ps}$",fontsize=18)
        plt.savefig("pos_vel_sort.png")

    if nd == 3:
        smooth_positions = np.zeros((len(positions)-2.0*n_avg,nd))
        times_vel = np.linspace(0.0,simulation_len,len(smooth_positions[skip_traj:,0]))
        times_pos = np.linspace(0.0,simulation_len,len(positions[:,0]))
        for i,val in enumerate(positions):
            if i > n_avg and i < len(positions)-2.0*n_avg:
                for ind in range(nd):
                    smooth_positions[i,ind] = np.mean(positions[i-n_avg:i+n_avg,ind])

        dy = step2ps
        print "Taking gradient of position data"
        velocities = np.zeros(np.shape(smooth_positions[skip_traj:,:]))
        for ind in range(nd):
            velocities[:,ind] = np.gradient(smooth_positions[skip_traj:,ind],dy)
        
        plt.clf()
        plt.figure()
        plt.subplot(212)
        plt.plot(times_vel[:],velocities[:,0],label="velocity")
        plt.ylabel(r"$\mathrm{Defect\ Velocity\ \AA/ps}$") 
        plt.xlabel(r"$\mathrm{Time\ ps}$") 

        plt.subplot(211)
        plt.plot(times_pos,positions[:,0],label="positions",alpha=0.2)
        plt.plot(times_vel,smooth_positions[skip_traj:,0],label="smooth positions")
        plt.ylabel(r"$\mathrm{Defect\ Position\ \AA}$") 
        plt.xlabel(r"$\mathrm{Time\ ps}$") 
        #plt.legend()
        plt.savefig(fn_velocities)

        plt.clf()
        plt.figure()
        plt.plot(times_pos[:20000], positions[10000:30000,0],alpha=0.2)
        plt.plot(times_vel[:20000], smooth_positions[10000:30000,0])
    #    plt.plot(times_vel[:20000], velocities[:20000,0])
        plt.savefig("test.png")

    print "Maximum absolute velocity = %14.5f Angstrom/ps" % (np.amax(abs(velocities)))
    print "Minimum absolute velocity = %14.5f Angstrom/ps" % (np.amin(abs(velocities)))

    return velocities

def test_normality(data,nbins=50.0):

    fn_norm_test = "norm_test.png"
   
    histvals1,xedges = np.histogram(data,nbins)

    print "Normalizing probability distribution"
    histvals = [(float(i)/float(sum(histvals1)))*100.0 for i in histvals1]
    xc = 0.5 * (xedges[:-1] + xedges[1:])

    plt.clf()
    plt.figure()
    plt.plot(xc,histvals)
    plt.xlabel(r"$\mathrm{Defect\ Velocity\ \AA/ps}$") 
    plt.ylabel(r"$\mathrm{Probability}$")
    plt.savefig(fn_norm_test)

    print "\n"
    
    return xc,histvals

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


def fit_func(func,xdata,ydata):
    
    len_simulation = 5000.0*11.0 # ps
    p0 = [0.0,0.1,2.0]
    popt,pcov = curve_fit(func,xdata,ydata,p0=p0)
    perr = np.sqrt(np.diag(pcov))

    print "Standard Deviation Error on Parameters:"
    for p in perr:
        print "%5.5E" % (p)

    plt.clf()
    plt.plot(xdata,func(xdata,popt[0],popt[1],popt[2]),label='Fitted Distribution')
    plt.plot(xdata,func(xdata,p0[0],p0[1],p0[2]),label='Guess Distribution')
    plt.plot(xdata,ydata,label='Actual Distribution')
    plt.legend()
    plt.savefig("fit_comparison.png")
    print "Optimized Fit Parameters: "
    print popt[0]
    print popt[1]
    print popt[2]

    return
    
def fit_compare(func1,func2,xdata,ydata):
    
    p0 = [0.0,0.1,2.0]
    popt1,pcov1 = curve_fit(func1,xdata,ydata,p0=p0)
    popt2,pcov2 = curve_fit(func2,xdata,ydata,p0=p0)
    perr1 = np.sqrt(np.diag(pcov1))
    perr2 = np.sqrt(np.diag(pcov2))

    print "Standard Deviation Error on First Fit Parameters:"
    for p in perr1:
        print "%5.5E" % (p)
    print "\n"

    print "Standard Deviation Error on Second Fit Parameters:"
    for p in perr2:
        print "%5.5E" % (p)
    print "\n"

    plt.clf()
    plt.subplot(211)
    plt.plot(xdata,func1(xdata,popt1[0],popt1[1],popt1[2]),label=r'$\mathrm{Normal\ fit}$')
    plt.plot(xdata,ydata)
    plt.ylabel(r"${\mathrm{Probability\ (\%)}}$",fontsize=16)
    plt.axis([-2.0,2.0,0.0,4.0])
    plt.legend(fontsize=16)

    plt.subplot(212)
    plt.plot(xdata,func2(xdata,popt2[0],popt2[1],popt2[2]),label=r'$\mathrm{Cauchy\ Fit}$')
    plt.plot(xdata,ydata)
    plt.xlabel(r"$\mathrm{Velocity\ \AA/ps}$",fontsize=16)
    plt.ylabel(r"${\mathrm{Probability\ (\%)}}$",fontsize=16)
    plt.axis([-2.0,2.0,0.0,4.0])
    plt.legend(fontsize=16)

    plt.savefig("fit_comparison.png")
    print "Optimized First Fit Parameters: "
    print popt1[0]
    print popt1[1]
    print popt1[2]
    print "\n"

    print "Optimized Second Fit Parameters: "
    print popt2[0]
    print popt2[1]
    print popt2[2]
    print "\n"

    return

def smooth_fft(
               velocities,
               ):

    fs2ps = 0.001
    stepsize = 50.0
    frate = 1.0/(stepsize*fs2ps)
    sp = np.fft.fft(velocities)
    freq = np.fft.fftfreq(velocities.shape[-1],frate)

    transfer_ps = 7.0
    transfer_step = transfer_ps*fs2ps*stepsize
    print transfer_step
    print frate
    print max(freq)/(stepsize*fs2ps)
    print min(abs(freq))

   # x = 0.00005
    x = 0.00002

    #x = 0.000005
    #x = 0.00001

    sp2 = np.copy(sp)
    freq2 = np.copy(freq)

    for i, val in enumerate(freq):
        if abs(val) > x:
            sp2[i] = 0.0

    smooth = np.real(np.fft.ifft(sp2))

    plt.clf()
    plt.figure()
    plt.plot(freq,sp,'b',alpha=0.3)
    plt.plot(freq2,sp2,'r',alpha=0.2)
    plt.axis([-2.0*x,2.0*x,0.0,0.2e8])
    plt.savefig("freq.png")

    end_sim = len(velocities)/frate
    x_axis = np.linspace(0.0,end_sim,len(smooth))
    x_axis2 = np.linspace(0.0,end_sim,len(velocities))

    plt.clf()
    plt.figure()
    plt.plot(x_axis2,velocities,'b',alpha=0.2)
    plt.plot(x_axis,smooth,'r',alpha=0.2)
    plt.savefig("irrev_freq.png")

    plt.clf()
    plt.figure()
    plt.plot(x_axis2[50:500],velocities[50:500],'b',alpha=0.3)
    plt.plot(x_axis[50:1000],smooth[50:1000],'g')
    plt.savefig("close_freq.png")

    return smooth
    
def test():
    
    nbins = 250.0

#    data = velocities(traj_len = 5000.0*11.0,
                      #ntraj = 1.0,
                      #nd = 1,
                      #axis = 2,
                      #fn_positions = "Ostar_all",
                      ##fn_positions = "Ostar",
                      #n_avg = 5,
                      #)
    #data = np.random.normal(0.0,0.1,1e6) #Tests normality Check
    #print "Average absolute velocity = %14.5f Angstrom/ps\n" % (np.mean(abs(data)))
    #xdata,ydata = test_normality(data,nbins)
    #fit_compare(gaussian,cauchy,xdata,ydata)

    fn_in = "Ostar10"
    nm2a = 10.0
    fs2ps = 0.001
    position = np.array(np.load(fn_in))[:,2]*nm2a
    position -= np.min(position)
    data = np.copy(position)
    smooth = smooth_fft(data)

    timestep = 50.0*fs2ps
    end = len(smooth)*timestep
    x = 5
    x2 = end

    t = np.linspace(0.0,end,len(smooth))
    velocity = np.gradient(smooth,timestep)

    acc = np.gradient(velocity,timestep)
    acc = np.real(acc)
    inds = []
    #for i, val in enumerate(acc):
    #for i, val in enumerate(velocity):
    #    if abs(val) < 1e-2:
    #        inds.append(i)
    #inds = np.array(inds)   
    positive_old = True
    for i, val in enumerate(velocity):
        if val > 0.0:
            positive_new = True
        else:
            positive_new = False
        if positive_new != positive_old:
            inds.append(i)
        positive_old = positive_new
    inds = np.array(inds)   

    dt = []
    dz = []
    v = []
    for n,ind in enumerate(inds):
        if n == 0:
            continue
        else:
            d = t[ind]-t[inds[n-1]]
            z = smooth[ind]-smooth[inds[n-1]]
            dt.append(d)
            dz.append(z)
            v.append(z/d)

    dt = np.array(dt)
    dz = np.array(dz)
    v = np.array(v)

    avg = np.average(abs(dz))
    std = np.std(abs(dz))
    dinds = []
    dinds2 = []
    dinds3 = []
    for p,delta in enumerate(dz):
        if abs(delta) > avg + 1.0*std:
            print delta
            dinds.append(inds[p])
            dinds2.append(inds[p+1])
            dinds3.append(inds[p-1])
    dinds = np.array(dinds) 
    dinds2 = np.array(dinds2)
    dinds3 = np.array(dinds3)
    

    plt.clf()
    plt.figure()
    plt.plot(t, data,alpha=0.3)
    plt.plot(t, smooth)
    #plt.plot(t[inds],smooth[inds],'ro')
    plt.plot(t[dinds],smooth[dinds],'go')
    plt.plot(t[dinds2],smooth[dinds2],'bo',alpha=0.5)
    #plt.plot(t[dinds3],smooth[dinds3],'ro',alpha=0.5)
    plt.axis([0.0,x2,0.0,max(data)+10])

    plt.savefig("TEST.png")

test()
    
