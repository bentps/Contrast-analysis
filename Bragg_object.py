#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 23 22:16:33 2017

@author: BenPS
"""
import matplotlib.pyplot as plt
from scipy.integrate import simps
from scipy.optimize import curve_fit
import numpy as np


#This functions are for use wtih the peak finding methods of the bragg class.
def my_cos(t):
    return np.cos(t)

def one_cos(x,shift,amp,width):
    y=np.piecewise(
        x, [(x>=0) & (x<shift),(x>=shift-width/2) & (x<(shift+width/2)),(x>=(shift+width))],
        [lambda x: 0, 
        lambda x: .5*amp*(-my_cos(2*np.pi*(x-shift+width/2)/width)+1), 
        lambda x: 0]
    )
    return y

#plt.plot(timebase,one_cos(timebase,1,.0001,.0002))

def cosfit(pair,shift,amp,width):
    p0=[shift,amp,width]
    timebase=pair[0]
    x=pair[1]
    fit = curve_fit(one_cos, timebase, x, p0=p0)
#    amplitude=fit[0][0]
#    error=np.sqrt(np.diag(fit[1]))[0]    
  
    return fit


#def my_cos(t):
#    return np.cos(t)
#
#def one_cos(x,shift,amp):
#    width=6.67294356543e-05
#    y=np.piecewise(
#        x, [(x>=0) & (x<shift),(x>=shift-width/2) & (x<(shift+width/2)),(x>=(shift+width))],
#        [lambda x: 0, 
#        lambda x: .5*amp*(-my_cos(2*np.pi*(x-shift+width/2)/width)+1), 
#        lambda x: 0]
#    )
#    return y
#
##plt.plot(timebase,one_cos(timebase,1,.0001,.0002))
#
#def cosfit(pair,shift,amp):
#    p0=[shift,amp]
#    timebase=pair[0]
#    x=pair[1]
#    fit = curve_fit(one_cos, timebase, x, p0=p0)
##    amplitude=fit[0][0]
##    error=np.sqrt(np.diag(fit[1]))[0]    
#  
#    return fit

class Bragg:
    def __init__(self, Phase):
        self.date=Phase.date
        self.N=Phase.N
        self.signals=Phase.Data.signals #dict: keys: twoT. Entries: tuples (timebase, signal)
        self.avg_signal=Phase.Data.avg_signal #dict: keys: twoT. Entries: tuple (timebase, avg_signal)
        self.time_list=Phase.Data.time_list #list of twoTs for this date, N.
        
        self.bragg_dict=Phase.Data.bragg_dict
        self.avg_bragg_dict=Phase.Data.avg_bragg_dict
        
        self.phase_dict=Phase.phase_dict#contains all phase info
        self.mean_phase_dict=Phase.mean_phase_dict# contains just phases for plotting, to be put into phaselist
        
        self.area_dict={}
        self.cor_phase_list=[]
        
        self.maxvolt_dict={}
        

        #Utility function for finding area under a curve
    def area(self, x):
        area = simps(x, dx=1)
        return area
        
    #Use this to find the bragg peak of interest. Once the appropriate 
    #first, last values are found, use them with bragg_area
    def plot_avg_bragg(self,twoT,first,last):
        triplet=self.avg_bragg_dict[twoT]
        print(self.area(triplet[1][first:last]))
        print(self.area(triplet[2][first:last]))        
        plt.plot(triplet[1][first:last],)
        plt.plot(triplet[2][first:last],)
#        np.savetxt('KDpulse.CSV', np.column_stack((triplet[0][first:last],triplet[1][first:last], triplet[2][first:last])), delimiter=",", fmt='%s')


    #This populates the area_dict for the selected bragg peak based on values
    #first, last
    def bragg_area(self, twoT, first, last):
        bragg_list=self.bragg_dict[twoT]
        if twoT in self.area_dict.keys():
            area_list=[]=self.area_dict[twoT]
        else:
            area_list=[]
            maxvolt1list=[]
            maxvolt2list=[]
        for i in bragg_list:
            area=self.area(i[0][1][first:last]*i[1][1][first:last])
#            area1=self.area(i[0][1][first:last])
#            area2=self.area(i[1][1][first:last])
            maxvolt1=max(i[0][1][first:last])
            maxvolt2=max(i[1][1][first:last])
#            self.area_list.append(area1*area2)
            area_list.append(area)
            maxvolt1list.append(maxvolt1)
            maxvolt2list.append(maxvolt2)
            
        self.area_dict[twoT]=area_list
        self.maxvolt_dict[twoT]=maxvolt1list,maxvolt2list
        
#        self.area_dict[twoT]=[i/2 for i in self.area_dict[twoT]]
        
        #quickly make plot of the area vs phase
    def area_vs_phase(self, twoT):        
        plt.plot(self.area_dict[twoT],self.phase_dict[twoT][0],'bo')
        plt.xlabel('Bragg pulse area^2 (arb units)')
        plt.ylabel('Phase (radians)')
        
        #This is for exporting the bragg data to send to Deep.
    def braggtracecsv(self, twoT):
        bragg_list=self.bragg_dict[twoT]
        timebase=bragg_list[0][0][0]
        braggtrace1=bragg_list[0][0][1]
        braggtrace2=bragg_list[0][1][1]
        braggtracecombined=bragg_list[0][0][1]*bragg_list[0][1][1]
        plt.plot(timebase,braggtrace1)
        plt.plot(timebase,braggtrace2)
        plt.plot(timebase,braggtracecombined)
#        np.savetxt('braggtracecombined.CSV', np.column_stack((timebase, braggtracecombined)), delimiter=",", fmt='%s')
        

        #functions for use in correct()
    def linear(self, x, A, B):
        return A*x + B
        
    def linearfit(self,x,y):
        fit = curve_fit(self.linear, x, y)
        return fit
        
    
    #Applies the simplist bragg correction and makes plots comparing
    #corrected to uncorrected
    def correct(self, twoT):
        area_list=self.area_dict[twoT]
        phase_list=self.phase_dict[twoT][0]
        fit=self.linearfit(area_list, phase_list)
        slope=fit[0][0]
        intercept=fit[0][1]
        print(slope)
        print(intercept)
        
#        base=list(np.linspace(min(area_list),max(area_list))) 
#        plt.plot(base,[self.linear(i,slope,intercept) for i in base])
#        plt.plot(area_list,phase_list,'ro')
#        plt.show
        
        combo=list(zip(area_list,phase_list))
        A0=np.mean(area_list)
        cor_phase_list=[i[1]-slope*(i[0]-A0) for i in combo]
        self.cor_phase_list=cor_phase_list
        plt.plot(cor_phase_list,'ro',label='corrected std ='+ str(np.std(cor_phase_list)))
        plt.plot(phase_list,'.',label='uncorrected std ='+ str(np.std(phase_list)))
        plt.ylabel('phase')
        plt.xlabel('shot number (N=1)')
        plt.legend()
        #plt.axis([0, 35, -.5, .5])

#        print(np.std(phase_list))
#        print(np.std(cor_phase_list))
#        plt.savefig('17-03-21_N1_correctedphases.pdf')     





        
    def ispeak(self, n, bragg):
        width=20
        mid=sum(bragg[n-width:n+width])
        left=sum(bragg[n-width-1:n+width-1])
        right=sum(bragg[n-width+1:n+width+1])
        if mid>=left and mid>=right:
            return True
        else:
            return False
        
    def peak_find(self,twoT):
#        x=self.avg_bragg_dict[twoT][1]
        x=self.avg_bragg_dict[twoT][1]*self.avg_bragg_dict[twoT][1]        
        peak_list=[]
#        plt.plot(x)
        for n in range(len(x)):
            if self.ispeak(n,x)==True:
                if n not in peak_list and n-1 not in peak_list and n+1 not in peak_list and x[n]>.01:
                    peak_list.append(n)
#        print(peak_list)
        return peak_list 
        
    def peak_space(self,twoT):
        x=self.avg_bragg_dict[twoT][1]*self.avg_bragg_dict[twoT][2]
        timebase=self.avg_bragg_dict[twoT][0]
        p=[timebase,x]
        peaklist=self.peak_find(twoT)
        timeguesslist=[timebase[i] for i in peaklist]
        timefitlist=[]
        timefiterrorlist=[]
        widthlist=[]
        for i in timeguesslist:
            fit=cosfit(p,i,.42,6.67e-05)
            timefit=fit[0][0]
            timefiterror=np.sqrt(np.diag(fit[1]))[0]
            width=fit[0][2]
            widthlist.append(width)
            timefitlist.append(timefit)
            timefiterrorlist.append(timefiterror)
        return timeguesslist, timefitlist, timefiterrorlist, widthlist
            
        
                

#    def peak_find(self,twoT):
#        x=list(self.avg_bragg_dict[twoT][1])
#        width=300
#        peak_list=[]
#        plt.plot(x)
#        for i in range(len(x)-width):
#            subset=x[i:i+width]
#            candidate=max(subset)
#            candidate_index=x.index(candidate)
#            if candidate_index not in peak_list and candidate>.2:
#                peak_list.append(candidate_index)
#        print(peak_list)
        


   