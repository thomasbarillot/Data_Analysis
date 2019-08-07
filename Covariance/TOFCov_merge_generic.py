# -*- coding: utf-8 -*-
"""
Created on Sun Jan 17 19:46:27 2016

@author: thomasbarillot
"""
import numpy as np
import scipy.io as sio
import glob

#import CovariancePlotting as cp
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--nos", help="enter number of samples",type=int)
parser.add_argument("--fpath", help="enter filepath")

args=parser.parse_args()


def CalculateCovMerge(arguments):
    
    fpath=arguments.fpath
    nsamples=arguments.nos
    
    #NOE=[]
    SX=np.zeros(2000,dtype=np.float64)
    SXI=np.zeros(2000,dtype=np.float64)
    SXX=np.zeros((2000,2000),dtype=np.float64)
    SI=np.zeros(1,dtype=np.float64)
    SI2=np.zeros(1,dtype=np.float64)
    
    #NOEfiles=glob.glob(fpath+'NOE_part*')
    SXfiles=glob.glob(fpath+'SX_part*')
    SXIfiles=glob.glob(fpath+'SXI_part*')
    SXXfiles=glob.glob(fpath+'SXX_part*')
    SIfiles=glob.glob(fpath+'SI_part*')
    SI2files=glob.glob(fpath+'SI2_part*')
    
    for i in range(0,len(SXfiles)):
        
        #NOE+=np.loadtxt(NOEfiles[i])
        SX+=np.loadtxt(SXfiles[i])
        SXI+=np.loadtxt(SXIfiles[i])
        SXX+=np.loadtxt(SXXfiles[i])
        SI+=np.loadtxt(SIfiles[i])
        SI2+=np.loadtxt(SI2files[i])
    
    CovElements={'SX':SX,'SXI':SXI,'SXX':SXX,'SI':SI,'SI2':SI2}
    
    sio.savemat(fpath+'CovElements.mat',CovElements)
    
CalculateCovMerge(args)

# Covariance
#    N=9.0*10**4
#    covX=(SXX-(np.outer(SX,SX)/N))/(N-1.0)
#    covXI=(SXI-(SX*SI)/N)/(N-1.0)
#    covII=(SI2-(SI**2)/N)/(N-1.0)

#    pcovX=((N-1.0)/(N-2.0))*(covX-(np.outer(covXI,covXI)/covII))
#covXI=(SXI/N-((SX*SI)/N**2)


#pathname=('%s%s/%s') %(folderpath,foldername,'pcovX.dat')
#np.savetxt(pathname,pcovX)
