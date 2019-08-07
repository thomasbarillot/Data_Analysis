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
parser.add_argument("--logpath", help="enter logfile pathname",type=str)

args=parser.parse_args()
sys.path.append('/Users/tbarillot/Documents/PostDoc_imperial_LabMac/007_repositories/SingleShotData/PeakFinding')
from PeakFinder import PeakFinder

from datetime import datetime
#import matplotlib.pyplot as plt
from mpi4py import MPI
comm=MPI.COMM_WORLD
rank=comm.Get_rank()
size=comm.Get_size()



def runPFFromSSData(arguments,logfile):
        
    noe=arguments.noe
    nsamples=arguments.nos
    path=arguments.fpath
    spath=arguments.spath
    
    ## Load template
    
    templatearray=np.load('/Users/tbarillot/Documents/PostDoc_imperial_LabMac/007_repositories/SingleShotData/PeakFinding/SingleShotTemplate.npz')
    template=templatearray['templateAmp']
    
    ##
    
    filelist=glob.glob(path+'specfile*.mat')
    #logfile.write("process scan %s rank %i\n" % (path,rank))
    #print len(filelist)
    
    count=0
    
    peaks_list=[]
    #amps_list=[]
    IR_list=[]
    failcount=0
    for i in range (0,len(filelist)):
            
        if i%(size)!=rank: continue # different ranks look at different file
        

        specstruct=sio.loadmat(path+'specfile_%s.mat' % str(i).zfill(3))
        #print specstruct.keys()

        specmattof=np.array(specstruct['specmat_ChA']).reshape((noe,nsamples))
        specmatIR=np.array(specstruct['specmat_ChB']).reshape((noe,nsamples))
        
        for j in range(0,noe):
            
            SStofiteration=-specmattof[j,:]
            #Reconstructtof=np.zeros((len(-specmattof[j,:])))
            convcriteria=1
            itercount=0
            
            peaks_list_tmp=[]
            #failcount=0
            try:
                
                
                while convcriteria != 0:

                    peaks,amps=PeakFinder(SStofiteration).StandardAnalysis()
                    peaks_list_tmp.append(peaks.astype(np.int16))
                    
                    
                    
                    for l in range(0,len(peaks)):
                        
                        
                        if peaks[l]>29 or peaks[l]<nsamples-120:
                     
                            tmptemplate=np.pad(template,(peaks[l]-30,nsamples-(peaks[l]+120)),'constant',constant_values=(0,0))
                            SStofiteration=SStofiteration-tmptemplate
                                
                            
                        elif peaks[l]<29:
                            try:
                                tmptemplate=np.pad(template[30-peaks[l]:],(0,nsamples-(peaks[l]+120)),'constant',constant_values=(0,0))
                                SStofiteration=SStofiteration-tmptemplate
                            except:
                                logfile.write('RANK%i Error substraction template down\n' % rank)
                        
                        elif peaks[l]>nsamples-120:
                            try:
                                tmptemplate=np.pad(template[:nsamples-peaks[l]+30],(peaks[l]-30,0),'constant',constant_values=(0,0))
                                SStofiteration=SStofiteration-tmptemplate
                            except:
                                logfile.write('RANK%i Error substraction template up\n' % rank)
            
                    baselineIR=np.median(specmatIR[j,nsamples-500:nsamples])
                    specmatIR[j,:]-=baselineIR
                    IR_list.append(specmatIR[j,300:335].sum().astype(np.int64))
                    
                    convcriteria=len(peaks)
                    itercount+=1
                
                    if itercount>200:
                        logfile.write('RANK%i WARNING: specfile %i shot %i reached the maximum of iterations without converging\n' % (rank,i,j))
                        break
                    
                peaks_list.append(np.concatenate(peaks_list_tmp))
            except:
                peaks_list.append(np.array([0],dtype=np.int16))
                failcount+=1
                #amps_list.append([0])
                IR_list.append(np.array([0],dtype=np.int16))
                logfile.write('RANK%i ERROR: Failed to process shot %i from specfile %i \n' % (rank,j,i))
                    
                
        
        count+=1
        #if rank==0:
        #    print 'Client 1:',count,'/',np.round(len(filelist)/8)
    
    if not os.path.exists(spath):
            os.makedirs(spath)
    
    peaks
    try:
        sio.savemat(spath+'SStofPeaks_part%s.mat' % str(rank).zfill(3),{'peaks':peaks_list,'IR':IR_list})
        logfile.write('\n Failed shot for rank %i : %i \n' % (rank,failcount))
        logfile.write('\n Client %i done \n' % rank)
    except:
        logfile.write('ERROR Saving: Client %i failed --> save locally \n' % rank)
        #sio.savemat()
        

dt=datetime
try:
    logfile=open(args.logpath,'a')
except:
    print  'ERROR'
logfile.write('\n launch peak finding processing  rank: %i, time %i%s%s_%s%s%s\n' % (rank,dt.today().year,str(dt.today().month).zfill(2),str(dt.today().day).zfill(2), \
                                        str(dt.today().hour).zfill(2), \
                                        str(dt.today().minute).zfill(2), \
                                        str(dt.today().second).zfill(2)))

runPFFromSSData(args,logfile)
MPI.Finalize()
logfile.close()
