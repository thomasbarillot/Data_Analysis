# -*- coding: utf-8 -*-
"""
Created on Tue Jul 25 09:16:23 2017

@author: Thomas Barillot
"""

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
import time

class CovMap():
    
    def __init__(self):
        
        
        ## Defines the CUDA kernel function
        self.mod = SourceModule("""
__global__ void XProductAndSumkernel(float *vectofX, float *vectofY, float *I, float *SX, float *SY, float *SXI, float *SYI, float *SXY)
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
                 SXI[Xidx_local]+=vectofX[Xidx_global]*I[k];
         }
         if(Xidx_local==0)
         {
                 SY[Yidx_local]+=vectofY[Yidx_global];
                 SYI[Yidx_local]+=vectofY[Yidx_global]*I[k];
         }
         
         SXY[EltIdx]+=vectofX[Xidx_global]*vectofY[Yidx_global];
     }
}
         
         
__global__ void XProductkernel(float *vectofX, float *vectofY, float *vectofXvectofY)
{

     int Xidx=threadIdx.x+blockIdx.x*blockDim.x;
     int Yidx=threadIdx.y+blockIdx.y*blockDim.y;
     int EltIdx=Xidx+Yidx*(gridDim.x*blockDim.x);
     
     vectofXvectofY[EltIdx]+=vectofX[Xidx]*vectofY[Yidx];

}
""")
    
        self.XProductAndSum = self.mod.get_function("XProductAndSumkernel")
        self.XProduct= self.mod.get_function("XProductkernel")
    
        ## Default values for covariance map parameters    
        self.noe=1000
        self.nsamples=8000
        self.nblocks=250
        self.bins=1
        self.fpath=''
        self.spath=''
        
        ## Allocate memory for the arrays
        self.NOE=np.zeros(1,dtype=np.float32)
        self.SX=np.zeros(self.nsamples,dtype=np.float32)
        self.SY=np.zeros(self.nsamples,dtype=np.float32)
        self.SXSY=np.zeros((self.nsamples,self.nsamples),dtype=np.float32)
        self.SXI=np.zeros(self.nsamples,dtype=np.float32)
        self.SYI=np.zeros(self.nsamples,dtype=np.float32)
        self.SXY=np.zeros((self.nsamples,self.nsamples),dtype=np.float32)
        self.SI=np.zeros(1,dtype=np.float32)
        self.SI2=np.zeros(1,dtype=np.float32)

        self.Cov=np.zeros((self.nsamples,self.nsamples),dtype=np.float32)
        self.pCov=np.zeros((self.nsamples,self.nsamples),dtype=np.float32)

# Calculate covariance elements in parallel
    def GetCovElements(self):
        
        self.NOE=[]
        self.SX=[]
        self.SY=[]
        self.SXSY=[]
        self.SXI=[]
        self.SYI=[]
        self.SXY=[]
        self.SI=[]
        self.SI2=[]
        self.Cov=[]
        self.pCov=[]
       
        self.NOE=np.zeros(1,dtype=np.float32)
        self.SX=np.zeros(self.nsamples,dtype=np.float32)
        self.SY=np.zeros(self.nsamples,dtype=np.float32)
        self.SXSY=np.zeros((self.nsamples,self.nsamples),dtype=np.float32)
        self.SXI=np.zeros(self.nsamples,dtype=np.float32)
        self.SYI=np.zeros(self.nsamples,dtype=np.float32)
        self.SXY=np.zeros((self.nsamples,self.nsamples),dtype=np.float32)
        self.SI=np.zeros(1,dtype=np.float32)
        self.SI2=np.zeros(1,dtype=np.float32)
        self.Cov=np.zeros((self.nsamples,self.nsamples),dtype=np.float32)
        self.pCov=np.zeros((self.nsamples,self.nsamples),dtype=np.float32)
        
        tmpSX=np.zeros(self.nsamples).astype(np.float32)
        tmpSY=np.zeros(self.nsamples).astype(np.float32)
        tmpSXI=np.zeros(self.nsamples).astype(np.float32)
        tmpSYI=np.zeros(self.nsamples).astype(np.float32)
        tmpSXY=np.zeros((self.nsamples,self.nsamples)).astype(np.float32)
    
        ## Compile the kernel
        
        fileslist=glob.glob(self.fpath+'specfile*.mat')
    
        nfiles=len(fileslist)    
    
        count=0
        print 'bp'
        for i in range(0,nfiles):
    
            specstruct=sio.loadmat(fileslist[i])
            #print self.nsamples,self.nsamples*self.noe
            tofmat=np.array(-specstruct['specmat_ChA'],order='C').astype(np.float32)
            #print tofmat.shape
            tmpI=np.reshape(tofmat,(1000,self.nsamples)).sum(1)
            
            
            self.XProductAndSum(cuda.In(tofmat),cuda.In(tofmat),cuda.In(tmpI),cuda.Out(tmpSX),cuda.Out(tmpSY),cuda.Out(tmpSXI),cuda.Out(tmpSYI),cuda.Out(tmpSXY),block=(32,32,1),grid=(self.nblocks,self.nblocks))
   
            self.SX+=tmpSX
            self.SY+=tmpSY
            self.SXY+=np.reshape(tmpSXY,(self.nsamples,self.nsamples))
        
            self.SI+=sum(tmpI)
            self.SI2+=sum(tmpI*tmpI)
            self.SXI+=tmpSXI
            self.SYI+=tmpSYI
        
            if i%(20)==0:
                print '%i'%i
        
            count+=1
        
        self.NOE=count*self.noe
        try:
            if not os.path.exists(self.spath):
                os.makedirs(self.spath)
        except:
            print 'Tried to create folder simultaneously for several nodes'

        self.XProduct(cuda.In(self.SX),cuda.In(self.SY),cuda.Out(self.SXSY),block=(32,32,1),grid=(self.nblocks,self.nblocks))
    
        CovElements={'SX':self.SX,'SY':self.SY,'SXSY':self.SXSY,'SXI':self.SXI,'SYI':self.SYI,'SXY':self.SXY,'SI':self.SI,'SI2':self.SI2,'NOE':self.NOE}
    
        np.savez(self.spath+'CovElements.npz',**CovElements)
        print 'merge done'


    def SetCovParameters(self,noe,nos,bins,fpath,spath):
        
        self.noe=noe
        
        self.nblocks=nos/32
        self.nsamples=self.nblocks*32
        self.bins=bins
        self.fpath=fpath
        self.spath=spath

    def GetCovMap(self,coeff_cov,pcov,coeff_pcov):
        
        self.Cov=(self.SXY-coeff_cov*self.SXSY/np.float32(self.NOE))/np.float32(self.NOE-1)
        
        if pcov==True:
            
            CovXI=((self.SXI-self.SX*self.SI/np.float32(self.NOE))/np.float32(self.NOE-1))
            CovYI=((self.SYI-self.SY*self.SI/np.float32(self.NOE))/np.float32(self.NOE-1))
            CovII=((self.SI2-self.SI*self.SI/np.float32(self.NOE))/np.float32(self.NOE-1))
            self.pCov=self.Cov-np.outer(CovXI,CovYI)/CovII
            
        
        







