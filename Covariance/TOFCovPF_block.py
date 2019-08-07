import numpy as np
import scipy.io as sio
import sys
import scipy.ndimage as im
import glob
import os

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--noe", help="enter number of events (shots) per files",type=int)
parser.add_argument("--nos", help="enter number of samples",type=int)
parser.add_argument("--binning", help="enter binning coeff",default=4.0,type=int)
parser.add_argument("--fpath", help="enter filepath")
parser.add_argument("--spath", help="enter savepath")
parser.add_argument("--xmin",default=0,type=int)
parser.add_argument("--xmax",default=8000,type=int)
args=parser.parse_args()

sys.path.append('/Users/tbarillot/Documents/PostDoc_imperial_LabMac/007_repositories/SingleShotData/PeakFinding/')
from PeakFinder import PeakFinder

from mpi4py import MPI
comm=MPI.COMM_WORLD
rank=comm.Get_rank()
size=comm.Get_size()



count=0;

#Fill the arrays with the sums

#for i in range(1,1):
    

def CalculateCov(arguments):
    
    noe=arguments.noe
    nsamples=arguments.nos
    fpath=arguments.fpath
    spath=arguments.spath
    xmin=arguments.xmin
    xmax=arguments.xmax
    binning=arguments.binning
    NOE=0
    
    SX=np.zeros((xmax-xmin)/binning,dtype=np.float64)
    SXI=np.zeros((xmax-xmin)/binning,dtype=np.float64)
    SXX=np.zeros(((xmax-xmin)/binning,(xmax-xmin)/binning),dtype=np.float64)
    SI=np.zeros(1,dtype=np.float64)
    SI2=np.zeros(1,dtype=np.float64)
    
    fileslist=glob.glob(fpath+'SStofPeaks_part*.mat')
    nfiles=len(fileslist)    
    
    count=0
    for i in range(0,nfiles):
    
        if i%(size)!=rank: continue # different ranks look at different file
        
        specstruct=sio.loadmat(fileslist[i])
        
        tofs=specstruct['peaks']
        #tofmat=np.reshape(tof,(noe,nsamples))[:,:6000]
        #tofmat*=-1
        
        for k in range(0,len(tofs[0])):
            
            peaksarray=tofs[0][k][0]
            tof,edges=np.histogram(peaksarray,bins=np.linspace(xmin,xmax,((xmax-xmin)/binning)+1))            
           
            tmpI=tof.sum()
            
            SXX+=np.outer(tof,tof)
            SXI+=tof*tmpI
            SI+=tmpI
            SI2+=tmpI**2
            SX+=tof
        
        count+=1
        
        if i%(size)==0:
            print '%i/%i' % (i,np.round(nfiles/4.0))
    
    NOE=noe*count
# Save the data in text files.        
    try:
        if not os.path.exists(spath):
                os.makedirs(spath)
    except:
        print 'Tried to create folder simultaneously for several nodes'
    
    
    savenameNOE='NOE_part%i.dat' % (rank)
    savenameSX='SX_part%i.dat' % (rank)
    savenameSXI='SXI_part%i.dat' % (rank)
    savenameSXX='SXX_part%i.dat' % (rank)
    savenameSI='SI_part%i.dat' % (rank)
    savenameSI2='SI2_part%i.dat' % (rank)

    #np.savetxt(spath+savenameNOE,NOE)
    np.savetxt(spath+savenameSX,SX)
    np.savetxt(spath+savenameSX,SX)
    np.savetxt(spath+savenameSXI,SXI)
    np.savetxt(spath+savenameSXX,SXX)
    np.savetxt(spath+savenameSI,SI)
    np.savetxt(spath+savenameSI2,SI2)


CalculateCov(args)
MPI.Finalize()







