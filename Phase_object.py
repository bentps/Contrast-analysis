# -*- coding: utf-8 -*-
"""
Created on Thu Apr 13 17:24:30 2017

@author: benps
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

import Data_object
import Phase_object

#Phase object takes Data object as an argument.

#changed phase dict to phase_dict dict 

class Phase:
    def __init__(self, Data):
        self.Data=Data
        self.date=Data.date
        self.N=Data.N
        self.signals=Data.signals #dict: keys: twoT. Entries: tuples (timebase, signal)
        self.avg_signal=Data.avg_signal #dict: keys: twoT. Entries: tuple (timebase, avg_signal)
        self.time_list=Data.time_list #list of twoTs for this date, N.
        
        self.phase_dict={}#contains all phase info
        self.mean_phase_dict={}# contains just phases for plotting, to be put into phaselist 
        self.std_dict={}
        self.std_error_dict={}
        
        self.phase_of_avg_dict={}
        self.phase_of_avg_amp={}
        self.phase_of_avg_offset={}
        
        self.ComplexAmpDict={}



#Run this before running lockfit to get envelope parameters
    def phaseampoffsetwidthfit(self, twoT, phase, amplitude, shift, offset, width):
        #this is a sine squared function
        def gauss(t, shift, width):
            return np.exp((-(t-shift)**2)/(2*width**2)) #width guess 0.00012
        def my_sin(t, phase, amplitude, offset):
            return ((np.sin(t*2*np.pi*4*omegarec + phase))**2)*amplitude+offset
        
        def envelope(t, phase, amplitude, shift, offset, width):
            return gauss(t, shift, width)*my_sin(t,phase,amplitude, offset)
        
        def sinefit(pair,phase, amplitude, shift, offset, width):
            p0=[phase, amplitude, shift, offset, width]
            timebase=pair[0]
            x=pair[1]
            fit = curve_fit(envelope, timebase, x, p0=p0)
            PhaseResult=fit[0][0]
            AmpResult=fit[0][1]
            ShiftResult=fit[0][2]
            OffsetResult=fit[0][3]
            WidthResult=fit[0][4]
            PhaseError=np.sqrt(np.diag(fit[1]))[0]
            AmpError=np.sqrt(np.diag(fit[1]))[1]
            ShiftError=np.sqrt(np.diag(fit[1]))[2]
            OffError=np.sqrt(np.diag(fit[1]))[3]
            WidthError=np.sqrt(np.diag(fit[1]))[4]
            
            data_fit = (timebase, envelope(timebase, *fit[0]))
          
            return PhaseResult, PhaseError, data_fit, AmpResult, OffsetResult, AmpError, OffError, WidthResult, WidthError, ShiftResult, ShiftError 

        avg=self.avg_signal[twoT]
        
        phase_of_avg=sinefit(avg, phase, amplitude, shift, offset, width)
        
        self.phase_of_avg_dict[twoT]=phase_of_avg
        
    def phaseampoffsetwidthfitall(self,phase, amplitude, shift, offset, width):
        for i in self.time_list:
            self.phaseampoffsetwidthfit(i,phase, amplitude, shift, offset, width) 
            
    def CompExp(self, t):
        return np.exp(-1j*16*np.pi*omegarec*t)

    def CompExp2(self, t):
        return np.exp(-2j*16*np.pi*omegarec*t)    
    
    
#I have two versions of lockfit. This version does a correction
#based on the evelope function. Therefore, it needs to have
#phaseampoffsetwidthfit run before it can run.    
    
#    def LockFit(self, twoT):
#        
#        #Gaussian function for use wtih getting envelope 
#        def gauss(t, shift, width, offset, amp):
#            return np.exp((-(t-shift)**2)/(2*width**2))*amp+offset #width guess 0.00012        
#        
#        ComplexAmpList=[]
#        PhaseList=[]
#        CorrectedComplexAmpList=[]
#        CorrectedPhaseList=[]        
#        for pair in self.signals[twoT]:
#            ExpList=[self.CompExp(t) for t in pair[0]]
#            ComplexAmp=np.dot(ExpList,pair[1])
#            ComplexAmpList.append(ComplexAmp)
#            
#        self.ComplexAmpDict[twoT]=ComplexAmpList
#
#        for amp in ComplexAmpList:
#            phase=(np.angle(amp)+np.pi)/2
#            PhaseList.append(phase)
#            
#        #Find average of all phases, for reliable phase guess during branch cut
#        PhaseOfMean=(np.angle(sum(ComplexAmpList))+np.pi)/2
##        print(PhaseOfMean)
#        #Center returned phase away from branch cut:
#        BranchPhaseList = [np.mod(phase-PhaseOfMean+np.pi/2,np.pi)+PhaseOfMean-np.pi/2 for phase in PhaseList]    
#        
#        #Get envelope info. This assumes that phaseampoffsetwidthfit has been run.
#        AvgSignalTuple=self.phase_of_avg_dict[twoT]
#        AvgShift=AvgSignalTuple[9]
#        AvgWidth=AvgSignalTuple[7]
#        AvgOffset=AvgSignalTuple[4]
#        AvgAmp=AvgSignalTuple[3]
#        
#        TimeBase=self.avg_signal[twoT][0]
#        
#        Envelope=[gauss(t,AvgShift,AvgWidth,AvgOffset,AvgAmp) for t in TimeBase]
#        
##        #Check that the evelope looks resonable, usually leave commented out
##        plt.plot(TimeBase,Envelope)
##        plt.plot(TimeBase,self.avg_signal[twoT][1])
#        
#        Exp0List=[1 for t in TimeBase]
#        Exp1List=[self.CompExp(t) for t in TimeBase]
#        Exp2List=[self.CompExp2(t) for t in TimeBase]
#        
#        E0=np.dot(Exp0List,Envelope)
#        E1=np.dot(Exp1List,Envelope)    
#        E2=np.dot(Exp2List,Envelope)
#        
#        Background=[0]*len(TimeBase) #Put shape of background noise if available
#        n1=np.dot(Exp1List,Background)
#        
#        Gamma=1 #Attenuation of transimpedence amp at 8*omegarec. Approximately 1.
#        
#
#
#        for amp in ComplexAmpList:
#            C0=amp/E0
#            Cold=C0
#            Cnew=(1/E0)*(amp-(2/Gamma)*np.abs(C0)*E1-np.conj(C0)*E2-n1)
#            while np.abs(Cnew-Cold)>1E-3:
#                Cold=Cnew
#                Cnew=(1/E0)*(amp-(2/Gamma)*np.abs(Cold)*E1-np.conj(Cold)*E2-n1)
##                print(Cold,Cnew)
#            
#            CorrectedComplexAmpList.append(Cnew)
#         
##        print(E0,E1,E2,np.angle(ComplexAmpList[0]),np.angle(CorrectedComplexAmpList[0]))    
#            
#        for amp in CorrectedComplexAmpList:
#            phase=(np.angle(amp)+np.pi)/2
#            CorrectedPhaseList.append(phase)
#            
#        #Find average of all phases, for reliable phase guess during branch cut
#        CorrectedPhaseOfMean=(np.angle(sum(CorrectedComplexAmpList))+np.pi)/2 
##        print(CorrectedPhaseOfMean)
#        CorrectedBranchPhaseList = [np.mod(phase-CorrectedPhaseOfMean+np.pi/2,np.pi)+CorrectedPhaseOfMean-np.pi/2 for phase in CorrectedPhaseList] 
#        
##        #use this if you want uncorrected phases
##        mean=np.mean(BranchPhaseList)
##        print(mean)
##        std=np.std(BranchPhaseList)
##        std_error=np.std(BranchPhaseList)/np.sqrt(len(BranchPhaseList))
##        self.phase_dict[twoT]=BranchPhaseList,[],[]  
##        self.mean_phase_dict[twoT]=mean
##        self.std_dict[twoT]=std                    
##        self.std_error_dict[twoT]=std_error
#
#        Weights=[np.abs(i) for i in CorrectedComplexAmpList]
##        print(Weights)
#        
#        #use this if you want corrected phases
#        mean=np.mean(CorrectedBranchPhaseList)
##        print(mean)
#        mean=np.average(CorrectedBranchPhaseList,weights=Weights)
#        std=np.std(CorrectedBranchPhaseList)
#        std_error=np.std(CorrectedBranchPhaseList)/np.sqrt(len(CorrectedBranchPhaseList))
#        self.phase_dict[twoT]=CorrectedBranchPhaseList,[],[]  
#        self.mean_phase_dict[twoT]=mean
#        self.std_dict[twoT]=std                    
#        self.std_error_dict[twoT]=std_error         

  
#This is the version of LockFit that can run without having an envelope
#function defined via phaseampoffsetwidthfit.
  
        
    def LockFit(self, twoT):
        ComplexAmpList=[]
        PhaseList=[]
        for pair in self.signals[twoT]:
            ExpList=[self.CompExp(t) for t in pair[0]]
            ComplexAmp=np.dot(ExpList,pair[1])
            ComplexAmpList.append(ComplexAmp)
            
        self.ComplexAmpDict[twoT]=ComplexAmpList

        for amp in ComplexAmpList:
            phase=(np.angle(amp)+np.pi)/2
            PhaseList.append(phase)
        
        PhaseOfMean=(np.angle(sum(ComplexAmpList))+np.pi)/2
        
        BranchPhaseList = [np.mod(phase-PhaseOfMean+np.pi/2,np.pi)+PhaseOfMean-np.pi/2 for phase in PhaseList]          

        mean=np.mean(BranchPhaseList)
        std=np.std(BranchPhaseList)
        std_error=np.std(BranchPhaseList)/np.sqrt(len(BranchPhaseList))
        self.phase_dict[twoT]=BranchPhaseList,[],[]  
        self.mean_phase_dict[twoT]=mean
        self.std_dict[twoT]=std                    
        self.std_error_dict[twoT]=std_error
        
        
  
            
                      
    def LockFitall(self):
        for i in self.time_list:
            self.LockFit(i)                       



            
            
            

#basic fitting method. Only phase is floated, no envelope, no variable amp, no offset.       
    def phasefit(self, twoT, phase, amp):
        #this is a sine squared function
        def my_sin(t, phase):
            amplitude=amp
            offset=0
            return ((np.sin(t*2*np.pi*4*omegarec + phase))**2)*amplitude + offset
        
        #fitting function. Check np.curve_fit documentation for info.
        #pair should be a tuple of (timebase, signal)
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
        
        #loop over signals at a given twoT
        signals=self.signals[twoT]
        phaselist=[]
        errorlist=[]
        fitlist=[]
        for i in signals:
            fit=sinefit(i, phase)
            phaselist.append(fit[0])
            errorlist.append(fit[1])
            fitlist.append(fit[2])
        mean=np.mean(phaselist)
        std=np.std(phaselist)
        std_error=np.std(phaselist)/np.sqrt(len(phaselist))
        self.phase_dict[twoT]=phaselist, errorlist, fitlist
        self.mean_phase_dict[twoT]=mean
        self.std_dict[twoT]=std                    
        self.std_error_dict[twoT]=std_error
        
        #get phase of avg at given twoT
        avg=self.avg_signal[twoT]
        avgfit=sinefit(avg, phase)
        phase_of_avg=avgfit[0]
        phase_of_avg_error=avgfit[1]
        phase_of_avg_datafit=avgfit[2]
        self.phase_of_avg_dict[twoT]=phase_of_avg, phase_of_avg_error, phase_of_avg_datafit
                                 
                                 
#Fit avg signal only      
    def avgphasefit(self, twoT, phase, amp):
        #this is a sine squared function
        def my_sin(t, phase):
            amplitude=amp
            offset=0
            return ((np.sin(t*2*np.pi*4*omegarec + phase))**2)*amplitude + offset
        
        #fitting function. Check np.curve_fit documentation for info.
        #pair should be a tuple of (timebase, signal)
        def sinefit(pair,phase):
            guess_phase = phase
            p0=[guess_phase]
            timebase=pair[0][0:2500]
            x=pair[1][0:2500]
            fit = curve_fit(my_sin, timebase, x, p0=p0)
            phase=fit[0][0]
            error=np.sqrt(np.diag(fit[1]))[0]    
            
            data_fit = (timebase, my_sin(timebase, *fit[0]))
          
            return phase, error, data_fit  
        
        
        #get phase of avg at given twoT
        avg=self.avg_signal[twoT]
        avgfit=sinefit(avg, phase)
        phase_of_avg=avgfit[0]
        phase_of_avg_error=avgfit[1]
        phase_of_avg_datafit=avgfit[2]
        self.phase_of_avg_dict[twoT]=phase_of_avg, phase_of_avg_error, phase_of_avg_datafit
        
    def phasefitall(self, phase, amp):
        for i in self.time_list:
            self.phasefit(i, phase, amp) 
            
    def avgphasefitall(self, phase, amp):
        for i in self.time_list:
            self.avgphasefit(i, phase, amp)             
  
    def phaseampoffsetfit(self, twoT, phase, amplitude, shift, offset):
        #this is a sine squared function
        def gauss(t, shift):
            return np.exp((-(t-shift)**2)/.00000003)
        def my_sin(t, phase, amplitude):
            return ((np.sin(t*2*np.pi*4*omegarec + phase))**2)*amplitude
        def envelope(t, phase, amplitude, shift, offset):
            return gauss(t, shift)*my_sin(t,phase,amplitude)+offset
        
        def sinefit(pair,phase, amplitude, shift, offset):
            guess_shift=shift
            guess_phase = phase
            guess_amp=amplitude
            guess_offset=offset
            p0=[guess_phase, guess_amp, guess_shift, guess_offset]
            timebase=pair[0][0:2500]
            x=pair[1][0:2500]
            fit = curve_fit(envelope, timebase, x, p0=p0, bounds=([-np.inf,-np.inf,-np.inf,-np.inf],[np.inf,np.inf,np.inf,np.inf]))
            phase=fit[0][0]
            amp=fit[0][1]
            shift=fit[0][2]
            offset=fit[0][3]
            error=np.sqrt(np.diag(fit[1]))[0]
            amperror=np.sqrt(np.diag(fit[1]))[1]
            offerror=np.sqrt(np.diag(fit[1]))[3]
            
            data_fit = (timebase, envelope(timebase, *fit[0]))
          
            return phase, error, data_fit, amp, offset, amperror, offerror  

        avg=self.avg_signal[twoT]
        
        phase_of_avg=sinefit(avg, phase, amplitude, shift, offset)
        
        self.phase_of_avg_dict[twoT]=phase_of_avg
                              
    def phaseampoffsetfitall(self,phase, amplitude, shift, offset):
        for i in self.time_list:
            self.phaseampoffsetfit(i,phase, amplitude, shift, offset)                              
        
        
    def phaseampoffsetrangefit(self, twoT, phase, amplitude, shift, offset,start,stop):
        #this is a sine squared function
        def gauss(t, shift):
            return np.exp((-(t-shift)**2)/.00000003)
        def my_sin(t, phase, amplitude):
            return ((np.sin(t*2*np.pi*4*omegarec + phase))**2)*amplitude
        def envelope(t, phase, amplitude, shift, offset):
            return gauss(t, shift)*my_sin(t,phase,amplitude)+offset
        
        def sinefit(pair,phase, amplitude, shift, offset):
            guess_shift=shift
            guess_phase = phase
            guess_amp=amplitude
            guess_offset=offset
            p0=[guess_phase, guess_amp, guess_shift, guess_offset]
            timebase=pair[0][start:stop]
            x=pair[1][start:stop]
            fit = curve_fit(envelope, timebase, x, p0=p0, bounds=([-np.inf,-np.inf,-np.inf,-np.inf],[np.inf,np.inf,np.inf,np.inf]))
            phase=fit[0][0]
            amp=fit[0][1]
            shift=fit[0][2]
            offset=fit[0][3]
            error=np.sqrt(np.diag(fit[1]))[0]
            amperror=np.sqrt(np.diag(fit[1]))[1]
            offerror=np.sqrt(np.diag(fit[1]))[3]
            
            data_fit = (timebase[start:stop], envelope(timebase[start:stop], *fit[0]))
          
            return phase, error, data_fit, amp, offset, amperror, offerror  

        avg=self.avg_signal[twoT]
        
        phase_of_avg=sinefit(avg, phase, amplitude, shift, offset)
        
        phase_of_avg_error=sinefit(avg, phase, amplitude, shift, offset)[1]
        phase_of_avg_datafit=sinefit(avg, phase, amplitude, shift, offset)[2]
        phase_of_avg_amp=sinefit(avg, phase, amplitude, shift, offset)[3]
        phase_of_avg_offset=sinefit(avg, phase, amplitude, shift, offset)[4]
        
        self.phase_of_avg_dict[twoT]=phase_of_avg  
        
        
#    def phaseampoffsetwidthfit(self, twoT, phase, amplitude, shift, offset, width):
#        #this is a sine squared function
#        def gauss(t, shift, width):
#            return np.exp((-(t-shift)**2)/(2*width**2)) #width guess 0.00012
#        def my_sin(t, phase, amplitude):
#            return ((np.sin(t*2*np.pi*4*omegarec + phase))**2)*amplitude
#        def envelope(t, phase, amplitude, shift, offset, width):
#            return gauss(t, shift, width)*my_sin(t,phase,amplitude)+offset
#        
#        def sinefit(pair,phase, amplitude, shift, offset, width):
#            p0=[phase, amplitude, shift, offset, width]
#            timebase=pair[0][0:2500]
#            x=pair[1][0:2500]
#            fit = curve_fit(envelope, timebase, x, p0=p0, bounds=([-np.inf,-np.inf,-np.inf,-np.inf],[np.inf,np.inf,np.inf,np.inf]))
#            PhaseResult=fit[0][0]
#            AmpResult=fit[0][1]
#            ShiftResult=fit[0][2]
#            OffsetResult=fit[0][3]
#            WidthResult=fit[0][4]
#            PhaseError=np.sqrt(np.diag(fit[1]))[0]
#            AmpError=np.sqrt(np.diag(fit[1]))[1]
#            OffError=np.sqrt(np.diag(fit[1]))[3]
#            WidthError=np.sqrt(np.diag(fit[1]))[4]
#            
#            data_fit = (timebase, envelope(timebase, *fit[0]))
#          
#            return PhaseResult, PhaseError, data_fit, AmpResult, OffsetResult, AmpError, OffError, WidthResult, WidthError 
#
#        avg=self.avg_signal[twoT]
#        
#        phase_of_avg=sinefit(avg, phase, amplitude, shift, offset, width)
#        
#        self.phase_of_avg_dict[twoT]=phase_of_avg

    def phaseampoffsetwidthfreqfit(self, twoT, phase, amplitude, shift, offset, width, freq):
        #this is a sine squared function
        def gauss(t, shift, width):
            return np.exp((-(t-shift)**2)/(2*width**2)) #width guess 0.00012
        def my_sin(t, phase, amplitude, freq):
            return ((np.sin(t*2*np.pi*4*freq + phase))**2)*amplitude
        def envelope(t, phase, amplitude, shift, offset, width, freq):
            return gauss(t, shift, width)*my_sin(t,phase,amplitude, freq)+offset
        
        def sinefit(pair,phase, amplitude, shift, offset, width, freq):
            p0=[phase, amplitude, shift, offset, width, freq]
            timebase=pair[0][0:2500]
            x=pair[1][0:2500]
            fit = curve_fit(envelope, timebase, x, p0=p0, bounds=([-np.inf,-np.inf,-np.inf,-np.inf],[np.inf,np.inf,np.inf,np.inf]))
            PhaseResult=fit[0][0]
            AmpResult=fit[0][1]
            ShiftResult=fit[0][2]
            OffsetResult=fit[0][3]
            WidthResult=fit[0][4]
            FreqResult=fit[0][5]
            PhaseError=np.sqrt(np.diag(fit[1]))[0]
            AmpError=np.sqrt(np.diag(fit[1]))[1]
            OffError=np.sqrt(np.diag(fit[1]))[3]
            WidthError=np.sqrt(np.diag(fit[1]))[4]
            FreqError=np.sqrt(np.diag(fit[1]))[5]
            
            data_fit = (timebase, envelope(timebase, *fit[0]))
          
            return PhaseResult, PhaseError, data_fit, AmpResult, OffsetResult, AmpError, OffError, WidthResult, WidthError, FreqResult, FreqError 

        avg=self.avg_signal[twoT]
        
        phase_of_avg=sinefit(avg, phase, amplitude, shift, offset, width, freq)
        
        self.phase_of_avg_dict[twoT]=phase_of_avg



    def plotvis(self):
        vis_list=[]
        time_list=self.time_list
        for i in self.phase_of_avg_offset:
            if self.phase_of_avg_offset[i]<0:
                self.phase_of_avg_offset[i]=0
        
        for i in time_list:
            amp=abs(self.phase_of_avg_amp[i]) # there is a hack on this line with abs(). To fix, need to restrict fit to positive amplitudes
            offset=self.phase_of_avg_offset[i]
            vis_list.append(amp/(amp+2*offset)) 
        plt.plot(time_list, vis_list,'ro', label='Visibility (unitless, range 0 to 1)')
        plt.legend(loc=0, fancybox=True, shadow=False, numpoints=1, frameon=True)
        plt.axis([0,14,0,3])
        plt.xlabel('2BigT (ms)')
        plt.ylabel('Arbitrary units')
          
    def showfits(self, twoT):
        signallist=self.signals[twoT]
        fitlist=self.phase_dict[twoT][2]
        for i in range(0,len(signallist)):
            data_fit=fitlist[i]
            signal=signallist[i]
            plt.plot(data_fit[0],data_fit[1], label='fit')
            plt.plot(signal[0],signal[1])
            plt.legend()
            plt.show()

          
            
    def show_fit_of_avg(self, twoT):
        signal=self.avg_signal[twoT]
        data_fit=self.phase_of_avg_dict[twoT][2]
        
#        amp=self.phase_of_avg_dict[.02][3]
#        offset=self.phase_of_avg_dict[.02][4]
#        vis=amp/(amp+2*offset)
        print(self.phase_of_avg_dict[twoT][0],self.phase_of_avg_dict[twoT][1])
        plt.plot(signal[0][0:2500],signal[1][0:2500],color='red',label='N='+str(self.N))
#        plt.plot(data_fit[0],data_fit[1], color='black', label='fit visibility = '+"%.3f" % vis)
        plt.plot(data_fit[0],data_fit[1], color='black')
#        plt.title('N='+str(self.N))
        plt.legend()
        plt.ylabel('PMT voltage',fontsize=14)
        plt.xlabel('seconds',fontsize=14)
#        plt.axis([0,.0005,0,1])
#        plt.savefig('N'+str(self.N)+'vis17-06-19.pdf')
        plt.show() 
        
        
    def show_fit_of_avg_all(self):
        for i in self.time_list:
            self.show_fit_of_avg(i)      

    def show_std_dev(self, twoT):
        phase_list=self.phase_dict[twoT][0]
        print(np.mean(phase_list), np.std(phase_list),self.std_dict[twoT], self.std_error_dict[twoT])
        plt.plot(phase_list,'ko')  
        plt.show()

#basic fitting method. Only phase is floated, no envelope, no variable amp, no offset.
#fitting to avgs of 10 signals       
    def clumpavgphasefit(self, twoT, phase, amp, Nsamp,Lsamp):
        #this is a sine squared function
        def my_sin(t, phase):
            amplitude=amp
            offset=0
            return ((np.sin(t*2*np.pi*4*omegarec + phase))**2)*amplitude + offset
        
        #fitting function. Check np.curve_fit documentation for info.
        #pair should be a tuple of (timebase, signal)
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
        
        #loop over signals at a given twoT
        signals=self.signals[twoT]
        Nsig=len(signals)
        phaselist=[]
        errorlist=[]
        fitlist=[]
        
        clumplist=[]
        avglist=[]
        
        
#loop over number of samples Nsamp
#for each sample, take Lsamp random signals and average them.         
        for i in range(Nsamp):
            sumdata=[0]*2500
            for j in range(Lsamp):
                rand=int(np.round((Nsig-1)*np.random.random()))
                randsig=signals[rand]
                randbase=signals[0][0]
                sumdata=sumdata+randsig[1]
                
            avglist.append((randbase,sumdata/Lsamp))    
                
        
        for i in avglist:
            fit=sinefit(i, phase)
            phaselist.append(fit[0])
            errorlist.append(fit[1])
            fitlist.append(fit[2])
        mean=np.mean(phaselist)
        std=np.std(phaselist)
        std_error=np.std(phaselist)/np.sqrt(len(phaselist))
        self.phase_dict[twoT]=phaselist, errorlist, fitlist
        self.mean_phase_dict[twoT]=mean
        self.std_dict[twoT]=std                    
        self.std_error_dict[twoT]=std_error
        
    def CompExpPhase(self, t, phase):
        return np.exp(-1j*(16*np.pi*omegarec*t+phase))
        
    def LockPhaseFit(self, twoT, phase):
        ComplexAmpList=[]
        PhaseList=[]
        for pair in self.signals[twoT]:
            ExpList=[self.CompExpPhase(t,phase) for t in pair[0]]
            ComplexAmp=np.dot(ExpList,pair[1])
            ComplexAmpList.append(ComplexAmp)
            
        self.ComplexAmpDict[twoT]=ComplexAmpList

        for amp in ComplexAmpList:
            Phase=(np.angle(amp)+np.pi)/2+phase/2
            PhaseList.append(Phase)
        
        PhaseOfMean=(np.angle(sum(ComplexAmpList))+np.pi)/2+phase/2
        
        BranchPhaseList = [np.mod(phi-PhaseOfMean+np.pi/2,np.pi)
                            +PhaseOfMean-np.pi/2 for phi in PhaseList]          
                       
#        mean=np.mean(BranchPhaseList)
        mean=np.mean(PhaseList)
        std=np.std(BranchPhaseList)
        std_error=np.std(BranchPhaseList)/np.sqrt(len(BranchPhaseList))
        return mean, std_error                   

    def LockPhase2pi(self, twoT):
        PhaseList=[]
        for phi in np.arange(0,2*np.pi,.1):
            Phase=self.LockPhaseFit(twoT, phi)[0]
            PhaseList.append(Phase)
        return PhaseList     