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
parser.add_argument("--fpath", help="enter filepath")
parser.add_argument("--spath", help="enter savepath")
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
    
    NOE=0
    SX=np.zeros(2000,dtype=np.float64)
    SXI=np.zeros(2000,dtype=np.float64)
    SXX=np.zeros((2000,2000),dtype=np.float64)
    SI=np.zeros(1,dtype=np.float64)
    SI2=np.zeros(1,dtype=np.float64)
    
    fileslist=glob.glob(fpath+'SStofPeaks_part*.mat')
    nfiles=len(fileslist)    
    
    count=0
    for i in range(0,nfiles):
    
        if i%(size)!=rank: continue # different ranks look at different file
        
        specstruct=sio.loadmat(fileslist[i])
        
        tofs=specstruct['peaks']
        
        for k in range(0,len(tofs[0,:])):
            

            peaksarray=tofs[0,k]
            if len(peaksarray)==0:
                continue
            tof,edges=np.histogram(np.concatenate(peaksarray),bins=np.linspace(0,8000,2001))            
           
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
    
    
    savenameNOE='NOE_part%i_%s.dat' % (rank,str(count).zfill(4))
    savenameSX='SX_part%i_%s.dat' % (rank,str(count).zfill(4))
    savenameSXI='SXI_part%i_%s.dat' % (rank,str(count).zfill(4))
    savenameSXX='SXX_part%i_%s.dat' % (rank,str(count).zfill(4))
    savenameSI='SI_part%i_%s.dat' % (rank,str(count).zfill(4))
    savenameSI2='SI2_part%i_%s.dat' % (rank,str(count).zfill(4))

    #np.savetxt(spath+savenameNOE,NOE)
    np.savetxt(spath+savenameSX,SX)
    np.savetxt(spath+savenameSX,SX)
    np.savetxt(spath+savenameSXI,SXI)
    np.savetxt(spath+savenameSXX,SXX)
    np.savetxt(spath+savenameSI,SI)
    np.savetxt(spath+savenameSI2,SI2)


CalculateCov(args)
MPI.Finalize()






