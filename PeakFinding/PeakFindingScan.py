# -*- coding: utf-8 -*-
"""
Created on Fri Jan 15 19:02:27 2016

@author: tbarillot
"""

import numpy as np
#import math
#import time
import scipy.io as sio
import scipy.ndimage as im
import sys
import glob
import os

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--noe", help="number of events(shots) per file",default=1000,type=int)
parser.add_argument("--nos", help="enter number of samples",type=int)
parser.add_argument("--fpath", help="enter filepath")
parser.add_argument("--spath", help="enter save filepath")

args=parser.parse_args()
sys.path.append('/Users/tbarillot/Documents/PostDoc_imperial_LabMac/007_repositories/SingleShotData/PeakFinding')
from PeakFinder import PeakFinder

from datetime import datetime
#import matplotlib.pyplot as plt
from mpi4py import MPI
comm=MPI.COMM_WORLD
rank=comm.Get_rank()
size=comm.Get_size()


def runAvgFromSSData(arguments):
        
    noe=arguments.noe
    nsamples=arguments.nos
    path=arguments.fpath
    spath=arguments.spath

        
    filelist=glob.glob(path+'specfile*.mat')
    print "process scan %s" % path
    #print len(filelist)
    
    count=0
    
    peaks_list=[]
    amps_list=[]
    IR_list=[]
    
    for i in range (0,len(filelist)):
            
        if i%(size)!=rank: continue # different ranks look at different file
        

        specstruct=sio.loadmat(path+'specfile_%s.mat' % str(i).zfill(3))
        #print specstruct.keys()

        specmattof=np.array(specstruct['specmat_ChA']).reshape((noe,nsamples))
        specmatIR=np.array(specstruct['specmat_ChB']).reshape((noe,nsamples))
        
        for j in range(0,noe):
            
            try:
                
                FTtof=np.fft.fft(-specmattof[j,:])
                if nsamples==8000:
                    FTtof[10:40]=im.median_filter(FTtof[10:40].real,3)
                    FTtof[7950:7990]=im.median_filter(FTtof[7950:7990].real,3)
                elif nsamples==5000:
                    FTtof[15:30]=im.median_filter(FTtof[15:30].real,3)
                    FTtof[4970:4985]=im.median_filter(FTtof[4970:4985].real,3)
                tmpcleanedtof=np.fft.ifft(FTtof)
                
                peaks,amps=PeakFinder(tmpcleanedtof.real).StandardAnalysis()
                peaks_list.append(peaks)
                amps_list.append(amps)
            
                baselineIR=np.median(specmatIR[j,nsamples-500:nsamples])
                specmatIR[j,:]-=baselineIR
                IR_list.append(specmatIR[j,300:335].sum())
            except:
                peaks_list.append([0])
                amps_list.append([0])
                IR_list.append([0])
                print 'Failed to process event %i from specfile %i' % (j,i)
        
        count+=1
        if rank==0:
            print 'Client 1:',count,'/',np.round(len(filelist)/8)
    
    if not os.path.exists(spath):
            os.makedirs(spath)
    
    sio.savemat(spath+'SStofPeaks_part%s.mat' % str(rank).zfill(3),{'peaks':peaks_list,'amps':amps_list,'IR':IR_list})
    print 'Client', rank, 'done'
        

dt=datetime
print 'launch peak finding processing  rank: %i, time %i%s%s_%s%s%s' % (rank,dt.today().year,str(dt.today().month).zfill(2),str(dt.today().day).zfill(2), \
                                        str(dt.today().hour).zfill(2), \
                                        str(dt.today().minute).zfill(2), \
                                        str(dt.today().second).zfill(2))

runAvgFromSSData(args)
MPI.Finalize()
