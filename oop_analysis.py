#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 26 00:13:07 2017

@author: BenPS
"""
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import pandas as pd
import os

#change the working directory to the data folder    
os.chdir('/Users/BenPS/Dropbox/Physics/Data') 

#this is the standard CSV import function, returning a numpy array of the oscilliscope voltage (signal) and time (timebase).
def csv(filename):
    the_data = pd.read_csv(filename, names=['a','b','c','d','e'])
    timebase=[i for i in the_data.c]
    timebase=np.array(timebase)
    signal=[i for i in the_data.d]
    signal=np.array(signal)
    return timebase, signal
    
#this is the avg function that should work for future data that hasn't been movingvanned
# date should be in quotes, in this format: '16-08-26'
#SorB should be 'signal' or 'bragg'.
#channel should be 'CH1' or 'CH2'. signal are generall all CH1, bragg data is on both channels.    
def avgofdate(date, SorB, channel, first, last, skiplist): 
    sumdata=[0]*2500
    for i in range(first,last):
        if i not in skiplist:
            if i<10:
                   filename=date+'/'+SorB+'/ALL000'+str(i)+'/F000'+str(i)+channel+'.CSV'
                   signal=csv(filename)[1]
                   sumdata=sumdata+signal 
            else:
                       
                filename=date+'/'+SorB+'/ALL00'+str(i)+'/F00'+str(i)+channel+'.CSV'
                signal=csv(filename)[1]
                sumdata=sumdata+signal  
    avg=sumdata/(last-first-len(skiplist))         
    return avg 

#function to round a number to the nearest multiple of pi, with prec digits.
def piround(x, prec=5, base=np.pi):
  return round(base * round(float(x)/base),prec)   

#plt.plot(avgofdate('17-03-03a','signal','CH1',20,40,[0]))


hbar=1.0545718E-34 #joule seconds
hbarerror=0.000000013E-34 #joule seconds

#myb=2.8883227E-22 #grams, from WolframAlpha

myb=2.8883227E-25 #kg, from WolframAlpha

#k=17992.007 #1/cm, from http://physics.nist.gov/PhysRefData/Handbook/Tables/ytterbiumtable5.htm
k=1799200.7 #1/m I'm pretty sure this is the vacuum value, but not positive.

omegarec=(2*np.pi*hbar*k**2)/(2*myb) #hz


#def linear(x, A, B):
#    return A*x + B
#    
#def linearfit(x,y,z):
#    fit = curve_fit(linear, x, y, sigma=z, absolute_sigma=True)
##    fit = curve_fit(linear, x, y)
#    return fit

class data:
    def __init__(self, date, N):
        self.date=date
        self.N=N
        self.time={}#sets of signals for a given time labled by time twoT.
        self.phase={}#contains all phase info
        self.avgsignal={}
        self.phase_of_avg_dict={}
        self.phase_of_avg_list=[]
        self.phase_of_avg_error_list=[]
        self.phase_of_avg_datafit={}
        self.phase_of_avg_amp={}
        self.mean_phase_dict={}# contains just phases for plotting, to be put into phaselist
        self.mean_phase_list=[]#phases for plotting
        self.std_error_dict={}
        self.std_error_list=[]
        self.time_list=[] #keeping track of
    
    def load(self, twoT, first, last, skiplist,letter):
        date=self.date
        self.time_list.append(twoT) #
        datapoints=[]
        sumdata=[0]*2500
        for i in range(first,last):
            if i not in skiplist:
                if i<10:
                    filename=date+letter+'/signal/ALL000'+str(i)+'/F000'+str(i)+'CH1.CSV'                        
                    signal=csv(filename)
                    datapoints.append(signal)
                    sumdata=sumdata+signal[1]                      
                else:
                    filename=date+letter+'/signal/ALL00'+str(i)+'/F00'+str(i)+'CH1.CSV'
                    signal=csv(filename)
                    datapoints.append(signal)
                    sumdata=sumdata+signal[1] 
        timebase=signal[0]            
        avg=sumdata/(last-first-len(skiplist)) 
        self.avgsignal[twoT]=(timebase, avg)
        self.time[twoT]=datapoints
    
    def phasefit(self, twoT, phase, amp):
        #this is a sine squared function
        def my_sin(t, phase):
            amplitude=amp
            offset=0
            return ((np.sin(t*2*np.pi*4*omegarec + phase))**2)*amplitude + offset
        def sinefit(pair,phase):
            guess_phase = phase
            p0=[guess_phase]
            timebase=pair[0]
            x=pair[1]
            fit = curve_fit(my_sin, timebase, x, p0=p0)
            phase=fit[0][0]
            error=np.sqrt(np.diag(fit[1]))[0]    
            
            data_fit = (timebase, my_sin(timebase, *fit[0]))
          
            return phase, error, data_fit  
        
        signals=self.time[twoT]
        phaselist=[]
        errorlist=[]
        fitlist=[]
        avg=self.avgsignal[twoT]
        phase_of_avg=sinefit(avg, phase)[0]
        phase_of_avg_error=sinefit(avg, phase)[1]
        phase_of_avg_datafit=sinefit(avg, phase)[2]
        self.phase_of_avg_dict[twoT]=phase_of_avg
        self.phase_of_avg_list.append(phase_of_avg)
        self.phase_of_avg_error_list.append(phase_of_avg_error)
        self.phase_of_avg_datafit[twoT]=phase_of_avg_datafit
        for i in signals:
            phaselist.append(sinefit(i, phase)[0])
            errorlist.append(sinefit(i, phase)[1])
            fitlist.append(sinefit(i, phase)[2])
        mean=np.mean(phaselist)
        std=np.std(phaselist)/len(phaselist)
        self.phase[twoT]=phaselist, errorlist, fitlist
        self.mean_phase_dict[twoT]=mean
        self.mean_phase_list.append(mean)
        self.std_error_dict[twoT]=std
        self.std_error_list.append(std)
        
    def phaseampfit(self, twoT, phase, amplitude, shift):
        #this is a sine squared function
        def gauss(t, shift):
            return np.exp((-(t-shift)**2)/.00000001)
        def my_sin(t, phase, amplitude):
            offset=0
            return ((np.sin(t*2*np.pi*4*omegarec + phase))**2)*amplitude + offset
        def envelope(t, phase, amplitude, shift):
            return gauss(t, shift)*my_sin(t,phase,amplitude)
        
        def sinefit(pair,phase, amplitude, shift):
            guess_shift=shift
            guess_phase = phase
            guess_amp=amplitude
            p0=[guess_phase, guess_amp, guess_shift]
            timebase=pair[0]
            x=pair[1]
            fit = curve_fit(envelope, timebase, x, p0=p0)
            phase=fit[0][0]
            amp=fit[0][1]
            shift=fit[0][1]
            error=np.sqrt(np.diag(fit[1]))[0]    
            
            data_fit = (timebase, envelope(timebase, *fit[0]))
          
            return phase, error, data_fit, amp  

        avg=self.avgsignal[twoT]
        phase_of_avg=sinefit(avg, phase, amplitude, shift)[0]
        phase_of_avg_error=sinefit(avg, phase, amplitude, shift)[1]
        phase_of_avg_datafit=sinefit(avg, phase, amplitude, shift)[2]
        phase_of_avg_amp=sinefit(avg, phase, amplitude, shift)[3]
        self.phase_of_avg_dict[twoT]=phase_of_avg
        self.phase_of_avg_list.append(phase_of_avg)
        self.phase_of_avg_error_list.append(phase_of_avg_error)
        self.phase_of_avg_datafit[twoT]=phase_of_avg_datafit    
        self.phase_of_avg_amp[twoT]=phase_of_avg_amp
        
    def phaseampfitall(self, phase, amp, shift):
        for i in self.time_list:
            self.phaseampfit(i, phase, amp, shift)


    def plotamp(self):
        amp_list=[]
        time_list=self.time_list
        for i in time_list:
            amp_list.append(self.phase_of_avg_amp[i])
        plt.plot(time_list, amp_list,'o')    
        
    def phasefitall(self, phase, amp):
        for i in self.time_list:
            self.phasefit(i, phase, amp)
        
    def showfits(self, twoT):
        signallist=self.time[twoT]
        fitlist=self.phase[twoT][2]
        for i in range(0,len(signallist)):
            data_fit=fitlist[i]
            signal=signallist[i]
            plt.plot(data_fit[0],data_fit[1], label='after fitting')
            plt.plot(signal[0],signal[1])
            plt.legend()
            plt.show()
                       
    def show_fit_of_avg(self, twoT):
        signal=self.avgsignal[twoT]
        data_fit=self.phase_of_avg_datafit[twoT]
        plt.plot(data_fit[0],data_fit[1], label='after fitting')
        plt.plot(signal[0],signal[1])
        plt.legend()
        plt.show() 
         
    #Function to adjust phases my multiples of pi to get them on the right slope. Apply to phasetimelists
    def phasepi(self):
            slope=93.321 #radians per millisecond
#            slope=93.32
    #        slope=121
            time_list=self.time_list
            mean_phase_list=self.mean_phase_list
            phase_of_avg_list=self.phase_of_avg_list
            N=self.N
            for n in range(0,len(mean_phase_list)):
                if n>0:
                    dt=time_list[n]-time_list[n-1]
                    dphistart=mean_phase_list[n]-mean_phase_list[n-1]
                    dphigoal=dt*slope*N**2
    #                dphigoal=dt*slope
                    addpi=piround(dphigoal-dphistart)
                    mean_phase_list[n]=mean_phase_list[n]+addpi      
            for n in range(0,len(phase_of_avg_list)):
                if n>0:
                    dt=time_list[n]-time_list[n-1]
                    dphistart=phase_of_avg_list[n]-phase_of_avg_list[n-1]
                    dphigoal=dt*slope*N**2
    #                dphigoal=dt*slope
                    addpi=piround(dphigoal-dphistart)
                    phase_of_avg_list[n]=phase_of_avg_list[n]+addpi                       
    def plot(self, whichphase):
        
        def linear(x, A, B):
            return A*x + B
        
        def linearfit(x,y,z):
            fit = curve_fit(linear, x, y, sigma=z, absolute_sigma=True)
        #    fit = curve_fit(linear, x, y)
            return fit

        time_list=self.time_list
        if whichphase==True:
            phase_list=self.mean_phase_list
            error_list=self.std_error_list
        if whichphase==False:
            phase_list=self.phase_of_avg_list
            error_list=self.phase_of_avg_error_list         
            
        N=self.N
        t = np.linspace(min(time_list), max(time_list), 256, endpoint=True)
#        t = np.linspace(.1, 2, 256, endpoint=True)
        fit=linearfit(time_list,phase_list,error_list)
        slope=fit[0][0]
        intercept=fit[0][1]
        error=np.sqrt(np.diag(fit[1]))[0] 
        
        omegaresult=slope*1000/(2*np.pi*4*N**2)
        errorresult=error*1000/(2*np.pi*4*N**2)
    
        plt.plot(t,linear(t,slope,intercept), label='N='+str(N))
#        plt.plot(t,linear(t,slope+error,intercept))
#        plt.plot(t,linear(t,slope-error,intercept))        
        plt.plot(time_list,phase_list,'r.')    
        plt.errorbar(time_list, phase_list, yerr=error_list, fmt='none', ecolor='black', elinewidth=1, capsize=4, capthick=1)
#        plt.axis([0,.25,-10,100])
        plt.legend(loc=0, fancybox=True, shadow=False, numpoints=1, frameon=True)
        plt.xlabel('2BigT (ms)')
        plt.ylabel('Phase (radians)')
#        print(error_list)
#        plt.show()
    
  

        print(['N='+str(N) , "%.3f" % slope, "%.3f" % error,  "%.3f" %omegaresult,  "%.3f" %errorresult])
        return ()                             
                             
        
#data1=data('17-03-03',10) #create data object instance
#
#data1.load(.040, 0, 10, [], '') #(twoT,firstfile,lastfile,skiplist, letter to append to date)
#data1.load(.042, 10, 20, [19], '')
#data1.load(.044, 20, 30, [], '')
#data1.load(.046, 30, 40, [], '')
#data1.load(.048, 40, 50, [], '')
#data1.load(.064, 50, 60, [], '')
#data1.load(.096, 60, 70, [], '')
#data1.load(.160, 70, 80, [], '')
#data1.load(.288, 80, 90, [], '')
#data1.load(.544, 90, 99, [90], '')
#data1.load(1.056, 0, 10, [7], 'a')
#data1.load(1.392, 10, 20, [10], 'a')
#data1.load(1.992, 20, 40, [], 'a')
#
#
##data1.phasefit(.04, 1, .3) #(twoT, phaseguess, sine amp)
##data1.phasefit(.042, 1, .3)
##data1.phasefit(.044, 1, .3)
##data1.phasefit(.046, 1, .3)
##data1.phasefit(.048, 1, .3)
#
#data1.phasefitall(0,.5)
#
#data1.phasepi()
#
#data1.plot(True)

#signal quality vs bigT vs atomnumber

#sigdata=data('17-03-06', 1)
#
#sigdata.load(2.0, 0, 5, [], '')
#sigdata.load(3.0, 5, 10, [], '')
#sigdata.load(4.0, 10, 15, [], '')
#sigdata.load(5.0, 15, 20, [], '')
#sigdata.load(6.0, 20, 25, [], '')
#sigdata.load(7.0, 25, 30, [], '')
#sigdata.load(8.0, 30, 35, [], '')
#sigdata.load(9.0, 35, 40, [], '')
#sigdata.load(1.0, 40, 45, [], '')
##sigdata.load(1.02, 45, 50, [], '')
##sigdata.load(1.04, 50, 55, [], '')
##sigdata.load(1.08, 55, 60, [], '')
##sigdata.load(1.16, 60, 65, [], '')
##sigdata.load(1.32, 65, 70, [], '')
#
#sigdata.phaseampfitall(0,.5, .0002)
#
#sigdata.plotamp()
#
#
#sigdata1=data('17-03-06', 1)
#
#sigdata1.load(1.0, 70, 75, [], '')
#sigdata1.load(3.0, 75, 80, [], '')
#sigdata1.load(5.0, 80, 85, [], '')
#sigdata1.load(7.0, 85, 90, [], '')
#sigdata1.load(9.0, 90, 99, [], '')
#
#sigdata1.phaseampfitall(0,.5, .0002)
#
##sigdata.show_fit_of_avg(7.0)
#
#sigdata1.plotamp()
#
#sigdata2=data('17-03-06', 1)
#
#sigdata2.load(1.0, 0, 5, [], 'a')
#sigdata2.load(2.0, 5, 10, [], 'a')
#sigdata2.load(3.0, 10, 15, [], 'a')
#sigdata2.load(4.0, 15, 20, [], 'a')
#sigdata2.load(5.0, 20, 25, [], 'a')
#
#sigdata2.phaseampfitall(0,.5, .0002)
#
#sigdata2.plotamp()
#
#plt.xlabel('2BigT (ms)')
#
#plt.ylabel('signal amplitude')
#
angdata1=data('17-03-08',1)

angdata1.load(1.0,15,20,[],'')
angdata1.load(2.0,20,25,[],'')
angdata1.load(3.0,25,30,[],'')
angdata1.load(4.0,30,35,[],'')
angdata1.load(5.0,35,40,[],'')

angdata1.phaseampfitall(0,.5, .0002)

angdata1.plotamp()

#angdata1.show_fit_of_avg(5.0)



#plt.savefig('100kblue_60korange_20kgreen_17-03-06.pdf')

#sigdata.phasefitall(1, .5)
##sigdata.phasefit(1.02, 0, .5)
#
#sigdata.phasepi()
#
#sigdata.plot(False)

#phaselist=data1.phase[.04][0]
#fitlist=data1.phase[.04][2]
#plt.plot(fitlist[0][0],fitlist[0][1])
#datapoints=data1.time[.04]
#plt.plot(datapoints[0][0],datapoints[0][1])

#data1.show_fit_of_avg(.04)
#data1.showfits(.04)
#plt.plot(data1.avgsignal[.04])

