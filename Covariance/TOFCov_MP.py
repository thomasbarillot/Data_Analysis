import numpy as np
import scipy.io as sio
import sys
import scipy.ndimage as im
import glob
import os

import multiprocessing as mp

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
def CalculateCov(arguments,rank):
    
    bins=arguments.bins
    noe=arguments.noe
    nsamples=arguments.nos
    
    fpath=arguments.fpath
    spath=arguments.spath
    
    NOE=np.zeros(1,dtype=np.float64)
    SX=np.zeros(nsamples/bins,dtype=np.float64)
    SXI=np.zeros(nsamples/bins,dtype=np.float64)
    SXX=np.zeros((nsamples/bins,nsamples/bins),dtype=np.float64)
    SI=np.zeros(1,dtype=np.float64)
    SI2=np.zeros(1,dtype=np.float64)
    
    fileslist=glob.glob(fpath+'specfile*.mat')
    
    # Get the current process id
    cores=arguments.cores
    
    nfiles=len(fileslist)    
    
    count=0
    for i in range(0,nfiles):
    
        if i%(cores)!=rank: continue # different ranks look at different file
        
        specstruct=sio.loadmat(fileslist[i])
        
        tofmat=np.int32(np.reshape(specstruct['specmat_ChA'],(noe,nsamples))[:,:nsamples])
        tofmat*=-1

        tofmat_binned=np.array(np.apply_along_axis(lambda a: np.histogram(np.arange(0,nsamples),bins=np.linspace(0,nsamples,nsamples/bins+1),weights=a)[0],1,tofmat))
       
       # tofmat_binned,ex,ey=np.histogram2d(np.arange(0,noe),np.arange(0,nsamples),bins=(np.arange(0,noe),np.linspace(0,nsamples,nsamples/bins+1)))
        
        SXX+=np.einsum('ij,ik->jk', tofmat_binned, tofmat_binned)
        SX+=tofmat_binned.sum(0)
        SXI+=np.einsum('ij,i->j', tofmat_binned, tofmat_binned.sum(1))
        SI+=sum(tofmat_binned.sum(1))
        SI2+=sum(tofmat_binned.sum(1)*tofmat_binned.sum(1))
        
        #for k in range(0,tofmat.shape[0]):
            
        #    tofvec=(tofmat[k,:])
            #PF=PeakFinder(tofvectmp)
            #PF.RemoveOffset(range(500,1000))
            #tofvec=np.array(PF.wf,dtype=np.float64)
            #tofvec[tofvec<200]=0
        #    tofvec_corrected,edges=np.histogram(np.arange(0,nsamples),bins=np.linspace(0,nsamples,nsamples/bins+1),weights=tofvec)
           
        #    tmpI=tofvec_corrected[:].sum()
            
        #    SXX=SXX+np.outer(tofvec_corrected,tofvec_corrected)
        #    SXI=SXI+tofvec_corrected*tmpI
        #    SI=SI+tmpI
        #    SI2=SI2+tmpI**2
        #    SX=SX+tofvec_corrected
        
        count+=1
        
        if i%(cores)==0:
            print '%i/%i' % (i/(cores),np.round(nfiles/cores))
    
    NOE+=noe*count
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

    np.savetxt(spath+savenameNOE,NOE)
    np.savetxt(spath+savenameSX,SX)
    np.savetxt(spath+savenameSX,SX)
    np.savetxt(spath+savenameSXI,SXI)
    np.savetxt(spath+savenameSXX,SXX)
    np.savetxt(spath+savenameSI,SI)
    np.savetxt(spath+savenameSI2,SI2)



# Merge the elements fron different processes
def CalculateCovMerge(arguments):
    
    print 'Start merging files...'
    bins=arguments.bins
    nsamples=arguments.nos
    fspath=arguments.spath
    cores=arguments.cores
    
    NOE=np.zeros(1,dtype=np.float64)
    SX=np.zeros(nsamples/bins,dtype=np.float64)
    SXI=np.zeros(nsamples/bins,dtype=np.float64)
    SXX=np.zeros((nsamples/bins,nsamples/bins),dtype=np.float64)
    SI=np.zeros(1,dtype=np.float64)
    SI2=np.zeros(1,dtype=np.float64)
    
    NOEfiles=glob.glob(fspath+'NOE_part*')
    SXfiles=glob.glob(fspath+'SX_part*')
    SXIfiles=glob.glob(fspath+'SXI_part*')
    SXXfiles=glob.glob(fspath+'SXX_part*')
    SIfiles=glob.glob(fspath+'SI_part*')
    SI2files=glob.glob(fspath+'SI2_part*')
    
    for i in range(0,cores):
        
        NOE+=np.loadtxt(NOEfiles[i])
        SX+=np.loadtxt(SXfiles[i])
        SXI+=np.loadtxt(SXIfiles[i])
        SXX+=np.loadtxt(SXXfiles[i])
        SI+=np.loadtxt(SIfiles[i])
        SI2+=np.loadtxt(SI2files[i])
        
        os.remove(NOEfiles[i])
        os.remove(SXfiles[i])
        os.remove(SXIfiles[i])
        os.remove(SXXfiles[i])
        os.remove(SIfiles[i])
        os.remove(SI2files[i])
        
    
    CovElements={'SX':SX,'SXI':SXI,'SXX':SXX,'SI':SI,'SI2':SI2,'NOE':NOE}
    
    np.savez(fspath+'CovElements.npz',**CovElements)
    print 'merge done'


#### Run the script ####

if __name__ == '__main__':

    ## Launch the file by file covariance compilation
    start=time.clock()
    print '%i cores start:'%argum.cores, start
    try:    
        processes=[mp.Process(target=CalculateCov,args=(argum,i)) for i in range(argum.cores)]
        print processes
    except:
        print 'ERROR def processes'


    for p in processes:
        try:
            p.start()
            print p
        except:
            print 'ERROR launch processes'
    for p in processes:
        print p.is_alive()
    for p in processes:
        p.join()
        
    ## Merge the data all together    
    intermstop=time.clock()
    print '%i cores merge start:'%argum.cores, intermstop
    try:
        CalculateCovMerge(argum)
    except:
        print 'ERROR merging files'
    
    stop=time.clock()
    print '%i cores stop:'%argum.cores, stop
        
    








