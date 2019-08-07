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
parser.add_argument("--x0",default=0,type=int)
parser.add_argument("--y0",default=0,type=int)

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
    x0=arguments.x0
    y0=arguments.y0
    binning=arguments.binning
    NOE=np.zeros(1,dtype=np.float64)
    SX=np.zeros(8000/binning,dtype=np.float64)
    SY=np.zeros(8000/binning,dtype=np.float64)
    SXI=np.zeros(8000/binning,dtype=np.float64)
    SYI=np.zeros(8000/binning,dtype=np.float64)
    SXY=np.zeros((8000/binning,8000/binning),dtype=np.float64)
    SI=np.zeros(1,dtype=np.float64)
    SI2=np.zeros(1,dtype=np.float64)
    
    fileslist=glob.glob(fpath+'specfile*.mat')
    nfiles=len(fileslist)    
    
    count=0
    for i in range(0,nfiles):
    
        if i%(size)!=rank: continue # different ranks look at different file
        
        specstruct=sio.loadmat(fileslist[i])
        
        tof=specstruct['specmat_ChA']
        tofmat=np.reshape(tof,(noe,nsamples))[:,:nsamples]
        tofmat*=-1
        
        for k in range(0,tofmat.shape[0]):
            
            tofvec=tofmat[k,:]
            #PF=PeakFinder(tofvectmp)
            #PF.RemoveOffset(range(500,1000))
            #tofvec=np.array(PF.wf,dtype=np.float64)
            #tofvec[tofvec<200]=0
            tofvec_corrected1,edges=np.histogram(np.arange(x0,x0+1000),bins=np.linspace(x0,x0+1000,(1000/binning)+1),weights=tofvec[x0:x0+1000])
            tofvec_corrected2,edges=np.histogram(np.arange(y0,y0+1000),bins=np.linspace(y0,y0+1000,(1000/binning)+1),weights=tofvec[y0:y0+1000])
            tmpI=tofvec.sum()
            
            SXY=SXY+np.outer(tofvec_corrected1,tofvec_corrected2)
            SXI=SXI+tofvec_corrected1*tmpI
            SYI=SYI+tofvec_corrected2*tmpI
            SI=SI+tmpI
            SI2=SI2+tmpI**2
            SX=SX+tofvec_corrected1
            SY=SY+tofvec_corrected2
        count+=1
        
        if i%(size)==0:
            print '%i/%i' % (i,np.round(nfiles/4.0))
    
    NOE+=noe*count
# Save the data in text files.        
    try:
        if not os.path.exists(spath):
                os.makedirs(spath)
    except:
        print 'Tried to create folder simultaneously for several nodes'
    
    
    savenameNOE='NOE_part%i.dat' % (rank)
    savenameSX='SX_part%i.dat' % (rank)
    savenameSY='SY_part%i.dat' % (rank)
    savenameSXI='SXI_part%i.dat' % (rank)
    savenameSYI='SYI_part%i.dat' % (rank)
    savenameSXY='SXY_part%i.dat' % (rank)
    savenameSI='SI_part%i.dat' % (rank)
    savenameSI2='SI2_part%i.dat' % (rank)

    np.savetxt(spath+savenameNOE,NOE)
    np.savetxt(spath+savenameSX,SX)
    np.savetxt(spath+savenameSY,SY)
    np.savetxt(spath+savenameSXI,SXI)
    np.savetxt(spath+savenameSYI,SYI)
    np.savetxt(spath+savenameSXY,SXY)
    np.savetxt(spath+savenameSI,SI)
    np.savetxt(spath+savenameSI2,SI2)

CalculateCov(args)
MPI.Finalize()







