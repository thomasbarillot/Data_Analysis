import numpy as np
import scipy.io
import scipy.stats.mstats 
#import matplotlib.pyplot as plt
#import IPython


import scipy.ndimage as im 
import math

convolv = scipy.ndimage.gaussian_filter
convolv1d = scipy.ndimage.gaussian_filter1d

class PeakFinder():
    def __init__(self,wf):
        self.wfraw=wf
        self.offset=0
        self.wf=wf
        self.rangelim=[0,len(self.wf)]

    def DefineRange(self,rangelim=[2000,-1]):
        if rangelim[1]==-1:
            rangelim[1]=len(self.wf)
        self.rangelim=rangelim
        self.wf=self.wf[rangelim[0]:rangelim[1]]

    def RemoveOffset(self,sfr=range(0,2000)):
        self.offset=np.median(self.wf[sfr]) #replaced mean by median
        self.wf=self.wf-self.offset
     
    def MedianFilter(self,points=3):
        self.wf = im.median_filter(self.wf,points)
    
    def NoiseThreshold(self,sfr=range(0,2000),factor=6):
        thr = factor*np.std(self.wf[sfr])
        self.wf = self.wf*(self.wf > thr)
        
    def GaussianFilter(self,width=12):
        self.wf = convolv1d(self.wf,width)
    
    def FindPeaks(self,extraconvolve=3):
        if(extraconvolve!=0):
            wfaux = convolv1d(self.wf,extraconvolve)
        else:
            wfaux = self.wf.copy()
        nzi = np.nonzero(wfaux[1:-2])[0] + 1  #nzi= non-zero indices
        rdiff = wfaux[1:] - wfaux[:-1]

        peaks = np.array([p for p in nzi if rdiff[p] < 0 and rdiff[p-1] > 0])+self.rangelim[0]
        return peaks
    
    def FindAllPeaksRefined(self,extraconvolve=3):
        peaks=self.FindPeaks(extraconvolve)       
        #peaks = np.array(map(self.QuadRefine, peaks))
        return peaks
    
    def FindPeaksNoRingingRefined(self,extend=400,factor=0.3,extraconvolve=3):
        #Ignores the peak that are within an "extend" factor from each peak and that are lower than peak*factor in high
        peaks=self.FindPeaks(extraconvolve)      
        validPeaks=np.zeros(len(peaks),dtype=np.bool)
        validPeaks[:]=True   #All peaks are valid at the beginning
        
        #For each peak in the list 
        for i in range(len(peaks)):
            threshold=self.wf[peaks[i]-self.rangelim[0]]*factor  #Define a threshold based on its amplitude
            for j in range(i+1,len(peaks)):     #And look for every subsequent peak
                if (peaks[j]<peaks[i]+extend):  #If it is within the extension of the previous peak
                    if(self.wf[peaks[j]-self.rangelim[0]]<threshold):  #And the peak high is smaller than the threshold
                        validPeaks[j]=False          #We mark the peak as not valid (potentially ringing)
                else:
                    break                       #As soon as we are outside the extension, we break the inner loop
                    
        
        #We save only the valid peaks        
        peaksNoRinging=[];
        for i in range(len(peaks)):
            if validPeaks[i]:
                peaksNoRinging.append(peaks[i])
                
           
        #Now we refine them
        peaksNoRinging = np.array(map(self.QuadRefinePlus, peaksNoRinging))
        return peaksNoRinging
               
    def QuadRefine(self,p):
        width=3
        pind=p-self.rangelim[0]
        minind=np.max([pind-width,0])
        maxind=np.min([len(self.wf),pind+width+1])
        x=np.arange(minind,maxind)+self.rangelim[0]
        y=self.wf[minind:maxind]
        coeffs=np.polyfit(x,y,2)

        pos=-coeffs[1] / (2*coeffs[0])
        
        if coeffs[0]<0 or pos>x[-1] or pos<x[0]:
            return p
        else:
            return pos
    
    def PeakAmplitudes(self,peaks,points=3):
        amps=[];
        for peak in peaks.round().astype(np.int32):
            amps.append(np.median(self.wfraw[peak-1:peak+1])-self.offset)

        amps=np.array(amps)

        return amps
        
        
    def StandardAnalysis(self):
        self.RemoveOffset(range(1000,2000));

        self.MedianFilter(2)
        self.GaussianFilter(1);
        self.NoiseThreshold(range(1000,2000),7);
        self.GaussianFilter(1);
        self.DefineRange([300,-1]);
        #peaks=self.FindPeaksNoRingingRefined(5,0.3)
        peaks=self.FindAllPeaksRefined(1)
        amps=self.PeakAmplitudes(peaks,3)
        
        inds=np.argsort(peaks)
        return peaks[inds],amps[inds]
    
    
        
def PeakFindTemplate(x,y,xtemplate,ytemplate,maxiterations=200):
 
    deadtime=xtemplate[-1] 

    goodpeaks=[]
    goodamps=[]
    signalaux=y.copy()
    iterations=maxiterations
    while(iterations>0):
        peaks,amps=PeakFinder(signalaux).StandardAnalysis()
        timepeaks=np.interp(peaks,range(len(x)),x)
        
        N=len(timepeaks)
        
        if N==0:
            break
        for i in range(N):
            if i>0 and (timepeaks[i]-timepeaks[i-1])<deadtime:
                continue
            signalaux=signalaux-np.interp(x,xtemplate+timepeaks[i],ytemplate,left=0, right=0)
            goodpeaks.append(timepeaks[i])
            goodamps.append(amps[i])
        iterations-=1

    if iterations==0:
        return None,None,None
        
    goodpeaks=np.array(goodpeaks)
    goodamps=np.array(goodamps)
        
    inds=np.argsort(goodpeaks)
    
    return goodpeaks[inds],goodamps[inds],signalaux

    
    
    
