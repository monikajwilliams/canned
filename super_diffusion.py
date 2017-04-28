#!/usr/bin/env python
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm
from numpy import linalg as LA
from scipy.optimize import curve_fit
import mdtraj as mdt

def t_alpha(x,scale,alpha):

    return scale*(x**alpha)

def t_alpha2(x,scale):

    return scale*(x**0.5)

def fit_func(func,
             xdata,
             ydata,
             p0=None,
             ):

    popt,pcov = curve_fit(func,xdata,ydata,p0=p0)
    popt2,pcov2 = curve_fit(t_alpha2,xdata,ydata)
    perr = np.sqrt(np.diag(pcov))

    #print "Standard Deviation Error on Parameters:"
    #for p in perr:
    #    print "%5.5E" % (p)

    plt.clf()
    plt.plot(xdata,func(xdata,popt[0],popt[1]),label='Fitted Scaling')
    plt.plot(xdata,t_alpha2(xdata,popt2[0]),label='Expected Scaling')
    # plt.plot(xdata,func(xdata,p0[0],p0[1]),label='Guess Distribution')
    plt.plot(xdata,ydata,label='Raw Data Scaling')
    plt.legend(loc=2)
    plt.savefig("fit_comparison.png")

    #print "Optimized Fit Parameters: "
    #print popt
    #print popt2


    return

def plot_rmsd(filename,
	   nd=1,
	   fn_in="msd_vals"	,
	   timestep=50,
	   fitrange=[0.0,90.0],
	   axis_range=0,
           cut = 0.0
           ):

        # => Loading data <= #
        if cut != 0.0:
            msd = np.load(fn_in)[:,:cut]
        else:
            msd = np.load(fn_in)
	times_steps = np.arange(len(msd[0]))

        # = > Unit Conversions <= #

	# Conversion factors
	nm22a2 = 100.0
	fs2ps = 0.001

	times = np.array(times_steps*timestep*fs2ps)
	msd_a = np.array(msd*nm22a2)
        rmsd_a = np.sqrt(msd_a)

        # setting up fit range
	nps = fitrange[0]
	nskip = int(nps/(fs2ps*timestep))
	eps = fitrange[1]
	eskip = int(eps/(fs2ps*timestep))

        # => Plotting <= #	

	plt.clf()
	plt.plot(np.log(times[1:]),np.log(rmsd_a[0,1:]),'ob', alpha=0.3,label='Data')


	if eskip-nskip >= len(times):
		m, b = np.polyfit(np.log(times[1:]), np.log(rmsd_a[0,1:]), 1)
		print "Skip range is longer than length of simulation!"
		D_value = False
	else:
		m, b = np.polyfit(np.log(times[nskip+1:eskip]), np.log(rmsd_a[0,nskip+1:eskip]), 1)
		D_value = True
                print "m = %f" % (m)

	plt.plot(np.log(times[nskip:eskip]),(np.log(times[nskip:eskip])*m)+b, alpha=1.0,c='r',label='Fitted Line')
        plt.xlabel(r'$\mathrm{log_{(t\ (ps))}}$',fontsize=20)
        plt.ylabel(r'$\mathrm{log_{(RMSD\ O^*\ (\AA))}}$',fontsize=20)
        plt.legend(loc=2)

	rmsd_max = max(rmsd_a[0])
	if axis_range != 0:
		plt.axis((axis_range[0],axis_range[1],axis_range[2],axis_range[3]))
	plt.savefig(filename) 

	if D_value == False:
		alpha = 0
		print "Error in simulation, diffusion coefficient is not valid!"
	else:
		alpha = m

        #fit_func(t_alpha,np.log(times[1:]),np.log(rmsd_a[0,1:]))
        #p0 = [b,m]
        #fit_func(t_alpha,times[1:],rmsd_a[0,1:],p0)

	#return alpha, b,np.log(rmsd_a[0,1:]),np.log(times[1:])
	return alpha

#def test():
#
#    alpha = plot_rmsd(filename='test.pdf')
#
#test()
