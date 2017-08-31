#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 23 22:10:37 2017

@author: BenPS
"""




import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

hbar=1.0545718E-34 #joule seconds
hbarerror=0.000000013E-34 #joule seconds
mybu=173.938867539 #unified atomic mass units, from ed myer
u=1.660539040E-27 #unified atomic mass units in kg
myb=mybu*u
vs=539386841661000 #hz
c=299792458 #m / s
k=vs/c #1/m
omegarec=(2*np.pi*hbar*k**2)/(2*myb) #hz


def linear(x, A, B):
        return A*x + B
    
def linearfit(x,y):
#    fit = curve_fit(linear, x, y, sigma=z, absolute_sigma=True)
    fit = curve_fit(linear, x, y)
    return fit

class Slope:
    def __init__(self, Phase):
        self.date=Phase.date
        self.N=Phase.N
        self.signals=Phase.Data.signals #dict: keys: twoT. Entries: tuples (timebase, signal)
        self.avg_signal=Phase.Data.avg_signal #dict: keys: twoT. Entries: tuple (timebase, avg_signal)
        self.time_list=Phase.Data.time_list #list of twoTs for this date, N.
        
        self.bragg_dict=Phase.Data.bragg_dict
        self.avg_bragg_dict=Phase.Data.avg_bragg_dict
        
        self.phase_dict=Phase.phase_dict#contains all phase info
        # contains just phases for plotting, to be put into phaselist     
        self.mean_phase_dict=Phase.mean_phase_dict   
        self.std_dict=Phase.std_dict
        self.std_error_dict=Phase.std_error_dict
        
        self.phase_of_avg_dict=Phase.phase_of_avg_dict
        self.phase_of_avg_amp=Phase.phase_of_avg_amp
        self.phase_of_avg_offset=Phase.phase_of_avg_offset
        
        self.mean_phase_list_pied=[]
        self.phase_of_avg_list_pied=[]
        self.residual_list=[]
        
        self.slope_fit=0
        self.phase_list=[]
        self.error_list=[]
        
#        self.omegalist=[]
#        self.omegaerror=[]
        
        self.guess_error=[]
        
    def piround(self, x, prec=5, base=np.pi):
        return round(base * round(float(x)/base),prec)       

    #Function to adjust phases my multiples of pi to get them on the right slope. Apply to phasetimelists
    def phasepi(self):
        slope=omegarec/(1000/(2*np.pi*4*1**2)) #radians per millisecond
#        slope=93
#        slope=121
        time_list=self.time_list
        mean_phase_list=[self.mean_phase_dict[twoT] for twoT in time_list]
        phase_of_avg_list=[self.phase_of_avg_dict[twoT][0] for twoT in time_list]
        N=self.N
        for n in range(0,len(mean_phase_list)):
            if n>0:
                dt=time_list[n]-time_list[n-1]
                dphistart=mean_phase_list[n]-mean_phase_list[n-1]
                dphigoal=dt*slope*N**2
#                dphigoal=dt*slope
                addpi=self.piround(dphigoal-dphistart)
                mean_phase_list[n]=mean_phase_list[n]+addpi      
        for n in range(0,len(phase_of_avg_list)):
            if n>0:
                dt=time_list[n]-time_list[n-1]
                dphistart=phase_of_avg_list[n]-phase_of_avg_list[n-1]
                dphigoal=dt*slope*N**2
#                dphigoal=dt*slope
                addpi=self.piround(dphigoal-dphistart)
                phase_of_avg_list[n]=phase_of_avg_list[n]+addpi  
                                 
        self.mean_phase_list_pied=mean_phase_list
        self.phase_of_avg_list_pied=phase_of_avg_list
        
        
    def avgphasepi(self):
        slope=omegarec/(1000/(2*np.pi*4*1**2)) #radians per millisecond
        time_list=self.time_list
        phase_of_avg_list=[self.phase_of_avg_dict[twoT][0] for twoT in time_list]
        N=self.N
     
        for n in range(0,len(phase_of_avg_list)):
            if n>0:
                dt=time_list[n]-time_list[n-1]
                dphistart=phase_of_avg_list[n]-phase_of_avg_list[n-1]
                dphigoal=dt*slope*N**2
#                dphigoal=dt*slope
                addpi=self.piround(dphigoal-dphistart)
                phase_of_avg_list[n]=phase_of_avg_list[n]+addpi  
                                 
#        self.mean_phase_list_pied=mean_phase_list
        self.phase_of_avg_list_pied=phase_of_avg_list        
                                    
    def phase_trial(self, whichphase):
            
        N=self.N 
#        slope=omegarec/(1000/(2*np.pi*4*1**2)) #radians per millisecond
        slope=3712/(1000/(2*np.pi*4*1**2)) #radians per millisecond
        
#        slope=93.321
#        slope=93.35
        time_list=self.time_list
        if whichphase=='mean_phase':
            phase_list=[self.mean_phase_dict[twoT] for twoT in time_list]
            error_list=[self.std_error_dict[twoT] for twoT in time_list]

            
            #base case
            base_list=phase_list[0:5]
            tbase_list=time_list[0:5]
            
            for n in range(0,len(base_list)):
                if n>0:
                    dt=tbase_list[n]-tbase_list[n-1]
                    dphistart=base_list[n]-base_list[n-1]
                    dphigoal=dt*slope*N**2
                    addpi=self.piround(dphigoal-dphistart)
                    base_list[n]=base_list[n]+addpi  
            
            for i in base_list: #apply base case to pied list
                self.mean_phase_list_pied.append(i)

                 
#            fit=linearfit(tbase_list,base_list)
#            base_slope=fit[0][0]
#            base_intercept=fit[0][1]
#            base_error=np.sqrt(np.diag(fit[1]))[0]  
#            t = np.linspace(min(tbase_list), max(tbase_list), 256, endpoint=True)    
#            plt.plot(t,linear(t,base_slope,base_intercept), label='N='+str(N))
#            plt.plot(tbase_list,self.mean_phase_list_pied,'r.')                  
            
            
            for phase in phase_list:
                index=phase_list.index(phase)
                if index > 4: #already dealt with first 5 phases in base case
                    sub_phase_list=self.mean_phase_list_pied[0:index]
                    sub_time_list=time_list[0:index]
                    
                    fit=linearfit(sub_time_list, sub_phase_list)
#                    slope=fit[0][0]
                    slope=93.321*N**2
#                    intercept=fit[0][1]
                    intercept=0
                    error=np.sqrt(np.diag(fit[1]))[0]
#                    print(slope)
                    
                    dt=time_list[index]-time_list[index-1]
#                    print(dt)
                    dphistart=-self.mean_phase_list_pied[index-1]+phase
#                    print(dphistart)                                   
                    dphigoal=dt*slope
#                    print(dphigoal)
                    addpi=self.piround(dphigoal-dphistart)
#                    print(addpi/np.pi)
                    self.mean_phase_list_pied.append(phase+addpi)
            
        if whichphase=='phase_of_avg':
            phase_list=[self.phase_of_avg_dict[twoT][0] for twoT in time_list]
            error_list=[self.phase_of_avg_dict[twoT][1] for twoT in time_list] 
            
            #base case
            base_list=phase_list[0:5]
            tbase_list=time_list[0:5]
            for n in range(0,len(base_list)):
                if n>0:
                    dt=tbase_list[n]-tbase_list[n-1]
                    dphistart=base_list[n]-base_list[n-1]
                    dphigoal=dt*slope*N**2
                    addpi=self.piround(dphigoal-dphistart)
                    self.guess_error.append(dphistart+addpi-dphigoal)                    
                    base_list[n]=base_list[n]+addpi

            
            for i in base_list: #apply base case to pied list
                self.phase_of_avg_list_pied.append(i)

            for phase in phase_list:
                index=phase_list.index(phase)
                if index > 4: #already dealt with first 5 phases in base case
                    sub_phase_list=self.phase_of_avg_list_pied[0:index]
                    sub_time_list=time_list[0:index]
                    
                    fit=linearfit(sub_time_list, sub_phase_list)
#                    slope=fit[0][0]
                    slope=93.321*N**2
                    intercept=0
#                    intercept=fit[0][1]
                    error=np.sqrt(np.diag(fit[1]))[0]
#                    print(slope)
                    
                    dt=time_list[index]-time_list[index-1]
#                    print(dt)
                    dphistart=-self.phase_of_avg_list_pied[index-1]+phase
#                    print(dphistart)                                   
                    dphigoal=dt*slope
#                    print(dphigoal)
                    addpi=self.piround(dphigoal-dphistart)
#                    print(addpi/np.pi)
                    self.phase_of_avg_list_pied.append(phase+addpi)                  
                    self.guess_error.append(dphistart+addpi-dphigoal)
#                
#        t = np.linspace(min(time_list), max(time_list), 256, endpoint=True)    
#        plt.plot(t,linear(t,base_slope,base_intercept), label='N='+str(N))
#        plt.plot(tbase_list,self.mean_phase_list_pied,'r.')  

                                 
    def plot(self, whichphase):
        
        def linear(x, A, B):
            return A*x + B
        
        def linearfit(x,y,z):
            fit = curve_fit(linear, x, y, sigma=z, absolute_sigma=True)
        #    fit = curve_fit(linear, x, y)
            return fit

        time_list=self.time_list
        if whichphase=='mean_phase':
            phase_list=self.mean_phase_list_pied
            error_list=[self.std_error_dict[twoT] for twoT in time_list]
        if whichphase=='phase_of_avg':
            phase_list=self.phase_of_avg_list_pied
            error_list=[self.phase_of_avg_dict[twoT][1] for twoT in time_list]     
            
        N=self.N

        phase_list=[i-phase_list[0] for i in phase_list]
        
        time_list=[i-time_list[0] for i in time_list]   

#        t = np.linspace(.1, 2, 256, endpoint=True)
        fit=linearfit(time_list,phase_list,error_list)
        slope=fit[0][0]
        intercept=fit[0][1]
        error=np.sqrt(np.diag(fit[1]))[0]
        
        omegaresult=slope*1000/(2*np.pi*4*N**2)#rename omega result
        errorresult=error*1000/(2*np.pi*4*N**2)
        
     
        fig = plt.figure(figsize=(4,7))

        ax1 = fig.add_subplot(211)
        ax2 = fig.add_subplot(212) 
        
        t = np.linspace(min(time_list), max(time_list), 1000, endpoint=True)    
        ax1.plot(t,linear(t,slope,0), label="%.3f" %omegaresult +"+/-"+  "%.3f" %errorresult)       
        ax1.plot(time_list,phase_list,'r.')    
        ax1.errorbar(time_list, phase_list, yerr=error_list, fmt='none', 
                     ecolor='black', elinewidth=1, capsize=4, capthick=1)
#        ax1.axis([0.0,.06,-10,600])
        ax1.legend(loc=0, fancybox=True, shadow=False, numpoints=1, frameon=True)
#        ax1.xlabel('2BigT (ms)')
#        ax1.ylabel('Phase (radians)')
        
        
        self.residual_list=[phase_list[i]-slope*time_list[i]-intercept for i in range(len(time_list))]
#        self.residual_list=[self.mean_phase_dict[twoT] for twoT in time_list]
        
        ax2.plot(time_list,[0]*len(time_list))
        ax2.plot(time_list, self.residual_list,'ro') 
        ax2.errorbar(time_list, self.residual_list, yerr=[1*i for i in error_list], 
                     fmt='none', ecolor='black', elinewidth=1, capsize=4, capthick=1)
#        ax2.axis([-.01,.03,-.1,.1])

        
        
        
#        self.residual_list=[phase_list[i]-slope*time_list[i]-intercept for i in range(len(time_list))]
##        self.residual_list=[self.mean_phase_dict[twoT] for twoT in time_list]
#        
#        plt.plot(time_list, self.residual_list,'ro') 
#        plt.errorbar(time_list, self.residual_list, yerr=error_list, 
#                     fmt='none', ecolor='black', elinewidth=1, capsize=4, capthick=1)
#

#        
#        t = np.linspace(min(time_list), max(time_list), 1000, endpoint=True)    
#        plt.plot(t,linear(t,slope,0), label='N='+str(N))
##        plt.plot(t,linear(t,slope+error,intercept))
##        plt.plot(t,linear(t,slope-error,intercept))        
##        plt.plot(time_list,phase_list,'r.')    
#        plt.errorbar(time_list, phase_list, yerr=error_list, fmt='none', 
#                     ecolor='black', elinewidth=1, capsize=4, capthick=1)
#        plt.axis([0.0,.06,-10,600])
#        plt.legend(loc=0, fancybox=True, shadow=False, numpoints=1, frameon=True)
#        plt.xlabel('2BigT (ms)')
#        plt.ylabel('Phase (radians)')
#        print(error_list)
#        plt.show()
    
#        omegalist.append(omegaresult)
#        omegaerror.append(errorresult)
#        
#        plt.errorbar(omegalist, self.N, yerr=omegaerror, fmt='none', 
#                     ecolor='black', elinewidth=1, capsize=4, capthick=1)
#        print(errorresult)

        print(['N='+str(N) , "%.5f" % slope, "%.5f" % error,
               "%.5f" %omegaresult,  "%.12f" %errorresult])
        return () 




    def getslope(self, whichphase):
        
        def linear(x, A, B):
            return A*x + B
        
        def linearfit(x,y,z):
            fit = curve_fit(linear, x, y, sigma=z, absolute_sigma=True)
        #    fit = curve_fit(linear, x, y)
            return fit

        time_list=self.time_list
        if whichphase=='mean_phase':
            phase_list=self.mean_phase_list_pied
            error_list=[self.std_error_dict[twoT] for twoT in time_list]
        if whichphase=='phase_of_avg':
            phase_list=self.phase_of_avg_list_pied
            error_list=[self.phase_of_avg_dict[twoT][1] for twoT in time_list]     
            
        N=self.N

#        t = np.linspace(.1, 2, 256, endpoint=True)
        fit=linearfit(time_list,phase_list,error_list)
        
        self.slope_fit=fit
        self.phase_list=phase_list
        self.error_list=error_list
        
        slope=fit[0][0]
        intercept=fit[0][1]
        error=np.sqrt(np.diag(fit[1]))[0]
        
        
        
#        #reset lists to start from zero
#        
#        phase_list=[i-phase_list[0] for i in phase_list]
#        
#        time_list=[i-time_list[0] for i in time_list]
#        
#        #plot residuals
##        self.residual_list=[phase_list[i]-slope*time_list[i] for i in range(len(time_list))]  
##        plt.plot(time_list, self.residual_list,'r.') 
##        plt.errorbar(time_list, self.residual_list, yerr=error_list, 
##                     fmt='none', ecolor='black', elinewidth=1, capsize=4, capthick=1)
#
#
#        
#        t = np.linspace(min(time_list), max(time_list)+1, 1000, endpoint=True)    
#        plt.plot(t,linear(t,slope,0), label='N='+str(N))
##        plt.plot(t,linear(t,slope+error,intercept))
##        plt.plot(t,linear(t,slope-error,intercept))        
#        plt.plot(time_list,phase_list,'r.')    
#        plt.errorbar(time_list, phase_list, yerr=error_list, fmt='none',
#                     ecolor='black', elinewidth=1, capsize=4, capthick=1)
#        plt.axis([0.0,.05,-10,600])
#        plt.legend(loc=0, fancybox=True, shadow=False, numpoints=1, frameon=True)
#        plt.xlabel('2BigT (ms)')
#        plt.ylabel('Phase (radians)')
##        print(error_list)
##        plt.show()
#    
#        omegaresult=slope*1000/(2*np.pi*4*N**2)
#        errorresult=error*1000/(2*np.pi*4*N**2)  
#
#        print(['N='+str(N) , "%.3f" % slope, "%.3f" % error,
#               "%.3f" %omegaresult,  "%.3f" %errorresult])
#        return () 