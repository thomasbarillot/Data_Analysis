import numpy as np
import scipy.io as sio
import sys
import scipy.ndimage as im
import glob
import os

import pycuda.autoinit
import pycuda.driver as cuda

from pycuda.compiler import SourceModule

from pycuda import gpuarray,tools


#import sys
#reload(sys)
#sys.setdefaultencoding('utf8')

import os

import time

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--cores", help="number of cores used",default=8,type=int)
parser.add_argument("--noe", help="enter number of events (shots) per files",type=int)
parser.add_argument("--nos", help="enter number of samples",type=int)
parser.add_argument("--bins", help="number of bins",type=int)
parser.add_argument("--fpath", help="enter filepath")
parser.add_argument("--spath", help="enter savepath")
argum=parser.parse_args()

#sys.path.append('/Users/tbarillot/Documents/PostDoc_imperial_LabMac/007_repositories/SingleShotData/PeakFinding/')

#from PeakFinder import PeakFinder

count=0;

# Calculate covariance elements in parallel
def CalculateCov(arguments,mod):
    
    mod = SourceModule("""
__global__ void Xproductkernel(float *vectofX, float *vectofY, float *SX, float *SY, float *SXY)
{
     for(int k=0;k<1000;k++){
         int Xidx_local=threadIdx.x+blockIdx.x*blockDim.x;
         int Xidx_global=threadIdx.x+blockIdx.x*blockDim.x+k*(gridDim.x*blockDim.x);
     
         int Yidx_local=threadIdx.y+blockIdx.y*blockDim.y;
         int Yidx_global=threadIdx.y+blockIdx.y*blockDim.y+k*(gridDim.y*blockDim.y);
     
         int EltIdx=Xidx_local+Yidx_local*(gridDim.x*blockDim.x);
         
         if(Yidx_local==0)
         {
                 SX[Xidx_local]+=vectofX[Xidx_global];
         }
         if(Xidx_local==0)
         {
                 SY[Yidx_local]+=vectofY[Yidx_global];
         }
         
         SXY[EltIdx]+=vectofX[Xidx_global]*vectofY[Yidx_global];
     }
}
""")
    
    
    
    bins=arguments.bins
    noe=arguments.noe
    nsamples=arguments.nos
    
    fpath=arguments.fpath
    spath=arguments.spath
    
    NOE=np.zeros(1,dtype=np.float32)
    SX=np.zeros(nsamples,dtype=np.float32)
    SY=np.zeros(nsamples,dtype=np.float32)
    
    SXI=np.zeros(nsamples,dtype=np.float32)
    SYI=np.zeros(nsamples,dtype=np.float32)
    SXY=np.zeros((nsamples,nsamples),dtype=np.float32)
    SI=np.zeros(1,dtype=np.float32)
    SI2=np.zeros(1,dtype=np.float32)
    
    tmpSX=np.zeros(nsamples).astype(np.float32)
    tmpSY=np.zeros(nsamples).astype(np.float32)
    tmpSXY=np.zeros((nsamples,nsamples)).astype(np.float32)
    
    ## Allocate memory on device for the matrix and sum arrays
    
   

    ## Compile the kernel
    
    Xproduct = mod.get_function("Xproductkernel")
    
    fileslist=glob.glob(fpath+'specfile*.mat')
    
    nfiles=len(fileslist)    
    
    count=0
    print 'bp'
    for i in range(0,nfiles):
    
        specstruct=sio.loadmat(fileslist[i])
        
        tofmat=np.array(-specstruct['specmat_ChA'],order='F').astype(np.float32)
        tmpI=tofmat.sum(1)

        Xproduct(cuda.In(tofmat),cuda.In(tofmat),cuda.Out(tmpSX),cuda.Out(tmpSY),cuda.Out(tmpSXY),block=(32,32,1),grid=(250,250))
        
        #cuda.memcpy_dtoh(tmpSX, SX_gpu)
        #cuda.memcpy_dtoh(tmpSY, SY_gpu)
        #cuda.memcpy_dtoh(tmpSXY, SXY_gpu)
        
        SX+=tmpSX
        SY+=tmpSY
        SXY+=np.reshape(tmpSXY,(nsamples,nsamples))
        
        SI+=sum(tmpI)
        SI2+=sum(tmpI*tmpI)
        SXI+=tmpSX*tmpI
        SYI+=tmpSY*tmpI
        
        if i%(20)==0:
            print '%i'%i
        
        count+=1
        
        
    try:
        if not os.path.exists(spath):
                os.makedirs(spath)
    except:
        print 'Tried to create folder simultaneously for several nodes'

    
    
    CovElements={'SX':SX,'SXI':SXI,'SYI':SYI,'SXY':SXY,'SI':SI,'SI2':SI2,'NOE':NOE}
    
    np.savez(spath+'CovElements.npz',**CovElements)
    print 'merge done'


#### Run the script ####

if __name__ == '__main__':

    
## kernel definition   


    
    ## Launch the covariance calculation function
    CalculateCov(argum,mod)
    








