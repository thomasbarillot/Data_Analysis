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
parser.add_argument("--x0",default=0,type=int)
parser.add_argument("--y0",default=8000,type=int)
parser.add_argument("--binning",default=4,type=int)
args=parser.parse_args()


def CalculateCovMerge(arguments):
    
    fpath=arguments.fpath
    nsamples=arguments.nos
    x0=arguments.x0
    y0=arguments.y0
    binning=arguments.binning

    NOE=np.zeros(1,dtype=np.float64)
    SX=np.zeros(1000/binning,dtype=np.float64)
    SY=np.zeros(1000/binning,dtype=np.float64)
    SXI=np.zeros(1000/binning,dtype=np.float64)
    SYI=np.zeros(1000/binning,dtype=np.float64)
    SXY=np.zeros((1000/binning,1000/binning),dtype=np.float64)
    SI=np.zeros(1,dtype=np.float64)
    SI2=np.zeros(1,dtype=np.float64)
    NOEfiles=glob.glob(fpath+'NOE_part*')
    SXfiles=glob.glob(fpath+'SX_part*')
    SYfiles=glob.glob(fpath+'SY_part*')
    SXIfiles=glob.glob(fpath+'SXI_part*')
    SYIfiles=glob.glob(fpath+'SYI_part*')
    SXYfiles=glob.glob(fpath+'SXY_part*')
    SIfiles=glob.glob(fpath+'SI_part*')
    SI2files=glob.glob(fpath+'SI2_part*')
    
    for i in range(0,len(SXfiles)):
        
        NOE+=np.loadtxt(NOEfiles[i])
        SX+=np.loadtxt(SXfiles[i])
        SY+=np.loadtxt(SYfiles[i])
        SXI+=np.loadtxt(SXIfiles[i])
        SYI+=np.loadtxt(SYIfiles[i])
        SXY+=np.loadtxt(SXYfiles[i])
        SI+=np.loadtxt(SIfiles[i])
        SI2+=np.loadtxt(SI2files[i])
    
    CovElements={'SX':SX,'SY':SY,'SXI':SXI,'SYI':SYI,'SXY':SXY,'SI':SI,'SI2':SI2,'NOE':NOE}
    
    sio.savemat(fpath+'CovElements_binning%i_x0_%i_y0_%i.mat'%(binning,x0,y0),CovElements)
    
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
