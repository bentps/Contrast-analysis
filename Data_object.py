# -*- coding: utf-8 -*-
"""
Created on Thu Mar 30 12:11:50 2017

@author: benps
"""

import pandas as pd
import numpy as np

#changed Time dict to signals dict

#hack to take subset of data. This is to make sure that he data consists of an
#integer number of oscillations.

#The contrast signals have 2500 voltage datapoints with .2us time step between data points.

#8/(.2E-6*8*omegarec)=1346.57   #peaks/(stepsize*freq) 
#2500-300-853=1347

#Use this if you want exactly 8 oscillations.
CutFromStart=300
CutFromEnd=853 

#Use this if you want the whole dataset. Note that usuing a non-integer number
#of oscillations leads to a systematic in the fitting.
#CutFromStart=0
#CutFromEnd=0 

class Data:
    def __init__(self, date, N):
        self.date=date
        self.N=N
        self.signals={}#sets of signals for a given time labled by time twoT.
        self.avg_signal={}
        self.time_list=[]
        self.bragg_dict={}
        self.avg_bragg_dict={}
#this is the standard CSV import function, returning a numpy array of the oscilliscope voltage (signal) and time (timebase).
    def csv(self, filename):
        the_data = pd.read_csv(filename, names=['a','b','c','d','e'])
        timebase=[i for i in the_data.c][CutFromStart:2500-CutFromEnd]
        timebase=np.array(timebase)
        signal=[i for i in the_data.d][CutFromStart:2500-CutFromEnd]
        signal=np.array(signal)
        return timebase, signal

    def load(self, twoT, first, last, skiplist, letter):
        date=self.date
        self.time_list.append(twoT) #
        if twoT in self.signals.keys():
            datapoints=self.signals[twoT]
        else:
            datapoints=[]
        if twoT in self.bragg_dict.keys():
            braggpoints=self.bragg_dict[twoT]
        else:
            braggpoints=[]
            
        if twoT in self.avg_signal.keys():
            sumdata=self.avg_signal[twoT][1]*(last-first-len(skiplist))
        else:
            sumdata=[0]*(2500-CutFromStart-CutFromEnd)
        if twoT in self.avg_bragg_dict.keys():
            sumbragg1=self.avg_bragg_dict[twoT][1]*(last-first-len(skiplist))
            sumbragg2=self.avg_bragg_dict[twoT][2]*(last-first-len(skiplist))
        else:
            sumbragg1=[0]*(2500-CutFromStart-CutFromEnd)
            sumbragg2=[0]*(2500-CutFromStart-CutFromEnd)         
        for i in range(first,last):
            if i not in skiplist:
                if i<10:
                    filename=date+letter+'/signal/ALL000'+str(i)+'/F000'+str(i)+'CH1.CSV'                        
                    signal_tuple=self.csv(filename)
                    datapoints.append(signal_tuple)
                    sumdata=sumdata+signal_tuple[1]  
                    
                    bragg1file=date+letter+'/bragg/ALL000{}/F000{}CH{}.CSV'.format(i,i,1)
                    bragg2file=date+letter+'/bragg/ALL000{}/F000{}CH{}.CSV'.format(i,i,2)
                    bragg1_tuple=self.csv(bragg1file)
                    bragg2_tuple=self.csv(bragg2file)
                    braggpoints.append((bragg1_tuple,bragg2_tuple))
                    sumbragg1=sumbragg1+bragg1_tuple[1]
                    sumbragg2=sumbragg2+bragg2_tuple[1]

                else:
                    filename=date+letter+'/signal/ALL00'+str(i)+'/F00'+str(i)+'CH1.CSV'                  
                    signal_tuple=self.csv(filename)
                    datapoints.append(signal_tuple)
                    sumdata=sumdata+signal_tuple[1] 
                    
                    bragg1file=date+letter+'/bragg/ALL00{}/F00{}CH{}.CSV'.format(i,i,1)
                    bragg2file=date+letter+'/bragg/ALL00{}/F00{}CH{}.CSV'.format(i,i,2)
                    bragg1_tuple=self.csv(bragg1file)
                    bragg2_tuple=self.csv(bragg2file)
                    braggpoints.append((bragg1_tuple, bragg2_tuple))
                    sumbragg1=sumbragg1+bragg1_tuple[1]
                    sumbragg2=sumbragg2+bragg2_tuple[1]
                    
        #avg signal stuff    
        timebase=signal_tuple[0]            
        avg=sumdata/(last-first-len(skiplist)) 
        self.avg_signal[twoT]=(timebase, avg)
        #signal stuff
        self.signals[twoT]=datapoints
        #avg bragg stuff
        braggbase=bragg1_tuple[0]    
        avg_bragg1=sumbragg1/(last-first-len(skiplist))
        avg_bragg2=sumbragg2/(last-first-len(skiplist))
        self.avg_bragg_dict[twoT]=(braggbase, avg_bragg1, avg_bragg2)
        #bragg stuff
        self.bragg_dict[twoT]=braggpoints
                    
                    
