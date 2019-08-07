import numpy as np
import scipy.io as sio
import sys
import scipy.ndimage as im
import glob
import os

import ctypes as c
import time
import sharedmem

## mp imports
import multiprocessing as mp
from multiprocessing import Pool,Lock

from scipy.stats import skewnorm
import scipy.optimize as opt

import time

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--cores", help="number of cores used",default=8,type=int)
parser.add_argument("--noe", help="enter number of events (shots) per files",type=int)
parser.add_argument("--nos", help="enter number of samples",type=int)
parser.add_argument("--bins", help="number of bins",type=int)
parser.add_argument("--fpath", help="enter filepath")
parser.add_argument("--fscan", help="enter filepath")
parser.add_argument("--spath",default='/media/thomasb/New Volume/Data/EPFL/NaCl_Auger/', help="enter savepath")
argum=parser.parse_args()

#sys.path.append('/Users/tbarillot/Documents/PostDoc_imperial_LabMac/007_repositories/SingleShotData/PeakFinding/')

#from PeakFinder import PeakFinder


# Calculate covariance elements in parallel
def CalculateCovWorker(arguments):
    
    startwork=time.time()
    
    bins=arguments[3]
    noe=arguments[1]
    nsamples=arguments[2]
    
    failed_evts=np.frombuffer(mp_failed_evts)
    NOE=np.frombuffer(mp_NOE)
    SX=np.frombuffer(mp_SX)
    SXI=np.frombuffer(mp_SXI)
    SXX=np.frombuffer(mp_SXX)
    SI=np.frombuffer(mp_SI)
    SI2=np.frombuffer(mp_SI2)
    
    
    # Get the current process id
    #cores=arguments.cores
    
    #nfiles=len(fileslist)    
    
    ## Fitting functions
    
    fitfunc = lambda p,x: p[0]*skewnorm.pdf(x,p[3],p[1],p[2])
    errfunc = lambda p,x,y: y-fitfunc(p,x)
    
    
        #if i%(cores)!=rank: continue # different ranks look at different file
        
    specstruct=sio.loadmat(arguments[0])
    
    tof=specstruct['specmat_ChA']
    I0raw = specstruct['specmat_ChB']
    tofmat=np.reshape(tof,(noe,nsamples))[:,:nsamples]
    tofmat*=-1
    I0raw_mat = -1*np.reshape(I0raw,(noe,nsamples))[:,:nsamples]
    
    blineI0=np.median(I0raw_mat[:,300:1000],axis=1)
    blinePES=np.median(tofmat[:,300:1000],axis=1)
    
    
    tofmat_corrected = (tofmat.T-blinePES).T[:,2300:3300]
    I0mat_corrected = (I0raw_mat.T-blineI0).T[:,1800:2500]
    
    count=0
    for k in range(0,tofmat_corrected.shape[0]):
        
        tofvec = np.int32(tofmat_corrected[k,:])
        tofvec[tofvec<500]=0
        I0vec = I0mat_corrected[k,:]
        
        p0 = [np.max(I0vec),np.argmax(I0vec),5,2]
        p1,ier = opt.leastsq(errfunc,p0,args=(np.arange(0,I0vec.shape[0],1),I0vec))
        #print ier
        if np.any(np.array([1,2,3,4])==ier)==False:
            #print p1
            lock.acquire()
            failed_evts+=1
            lock.release()
            continue
        
        fit=p1[0]*skewnorm.pdf(np.arange(0,I0vec.shape[0],1),p1[3],p1[1],p1[2])
        mask=np.zeros((I0vec.shape[0])).astype(np.bool)
        mask[fit>0.1*np.max(fit)]=1
        maskidx=np.where(mask==True)[0]
        
        try:
            SNR1=np.mean(I0vec[mask])/np.std(I0vec[maskidx[0]-200:maskidx[0]])
            SNR2=np.mean(I0vec[mask])/np.std(I0vec[maskidx[-1]:maskidx[-1]+200])
            if SNR1<6:
                lock.acquire()
                failed_evts+=1
                lock.release()
                continue
            if SNR2<6:
                lock.acquire()
                failed_evts+=1
                lock.release()
                continue
        except:
            lock.acquire()
            failed_evts+=1
            lock.release()
            
        I0val = skewnorm.std(p1[3],p1[1],p1[2])*p1[0]
        #I0val=np.sum(I0vec[mask])
        #PF=PeakFinder(tofvectmp)
        #PF.RemoveOffset(range(500,1000))
        #tofvec=np.array(PF.wf,dtype=np.float64)
        #tofvec[tofvec<200]=0
        tofvec_binned,edges=np.histogram(np.arange(0,nsamples/4),bins=np.linspace(0,nsamples/4,(nsamples/bins)/4+1),weights=tofvec)
        
        #tofvec_binned,edges=np.histogram(np.arange(0,nsamples),bins=np.linspace(0,nsamples,nsamples/bins+1),weights=tofvec)
        
        #tmpI=tofvec_corrected[:].sum()
        lock.acquire()
        SXX += np.concatenate(np.outer(tofvec_binned,tofvec_binned))
        SXI += tofvec_binned*I0val
        SI += I0val
        SI2 += I0val**2
        SX += tofvec_binned
        NOE += 1
        lock.release()
        
    print arguments[0]+'....Done in : ',time.time()-startwork

# Save the data in text files.        
    



# Merge the elements fron different processes
#def CalculateCovMerge(arguments):
#    
#    print 'Start merging files...'
#    bins=arguments.bins
#    nsamples=arguments.nos
#    fspath=arguments.spath
#    cores=arguments.cores
#    
#    NOE=np.zeros(1,dtype=np.float64)
#    SX=np.zeros(nsamples/bins,dtype=np.float64)
#    SXI=np.zeros(nsamples/bins,dtype=np.float64)
#    SXX=np.zeros((nsamples/bins,nsamples/bins),dtype=np.float64)
#    SI=np.zeros(1,dtype=np.float64)
#    SI2=np.zeros(1,dtype=np.float64)
#    
#    NOEfiles=glob.glob(fspath+'NOE_part*')
#    SXfiles=glob.glob(fspath+'SX_part*')
#    SXIfiles=glob.glob(fspath+'SXI_part*')
#    SXXfiles=glob.glob(fspath+'SXX_part*')
#    SIfiles=glob.glob(fspath+'SI_part*')
#    SI2files=glob.glob(fspath+'SI2_part*')
#    
#    for i in range(0,cores):
#        
#        NOE+=np.loadtxt(NOEfiles[i])
#        SX+=np.loadtxt(SXfiles[i])
#        SXI+=np.loadtxt(SXIfiles[i])
#        SXX+=np.loadtxt(SXXfiles[i])
#        SI+=np.loadtxt(SIfiles[i])
#        SI2+=np.loadtxt(SI2files[i])
#        
#        os.remove(NOEfiles[i])
#        os.remove(SXfiles[i])
#        os.remove(SXIfiles[i])
#        os.remove(SXXfiles[i])
#        os.remove(SIfiles[i])
#        os.remove(SI2files[i])
#        
#    
#    CovElements={'SX':SX,'SXI':SXI,'SXX':SXX,'SI':SI,'SI2':SI2,'NOE':NOE}
#    
#    np.savez(fspath+'CovElements.npz',**CovElements)
#    print 'merge done'


#### Run the script ####

if __name__ == '__main__':

    ## Launch the file by file covariance compilation
    
    
    start=time.time()
    
    Speclength = np.int(argum.nos/argum.bins)/4
    fpath=argum.fpath
    spath=argum.spath
    fscan=argum.fscan
    
    try:
        if not os.path.exists(spath):
                os.makedirs(spath)
    except:
        print 'Error making directory'
    
    lock=Lock()
    
    failed_evts = np.zeros(1,dtype=np.float64)
    NOE_tot = np.zeros(1,dtype=np.float64)
    SX_tot = np.zeros(Speclength,dtype=np.float64)
    SXI_tot = np.zeros(Speclength,dtype=np.float64)
    SXX_tot = np.zeros((Speclength,Speclength),dtype=np.float64)
    SI_tot = np.zeros(1,dtype=np.float64)
    SI2_tot = np.zeros(1,dtype=np.float64)
    
    Cov_Sliced= np.zeros((Speclength,Speclength),dtype=np.float64)
    pCov_Sliced= np.zeros((Speclength,Speclength),dtype=np.float64)
    
    mp_failed_evts = mp.Array(c.c_long,1,lock=False)
    mp_NOE = mp.Array(c.c_long,1,lock=False)
    mp_SX = mp.Array(c.c_long,Speclength,lock=False)
    mp_SXI = mp.Array(c.c_long,Speclength,lock=False)
    mp_SXX = mp.Array(c.c_long,Speclength*Speclength,lock=False)
    mp_SI = mp.Array(c.c_long,1,lock=False)
    mp_SI2 = mp.Array(c.c_long,1,lock=False)
   
    
    nfilesstop=180
    fileslist=glob.glob(fpath+'specfile*.mat')[:nfilesstop]
    slicestep=40
    slices=np.arange(0,len(fileslist),slicestep)    

    #pool=Pool(processes=7)
    
    for nf in slices:
        
        pool=Pool(processes=7)
        
        arg=[(f,argum.noe,argum.nos,argum.bins) for f in fileslist[nf:nf+slicestep]]
        res=pool.map(CalculateCovWorker,arg)
    #print res
        pool.close()
        pool.join()
        
    	
    	#ptr of the shared arrays 
        failed_evts=np.frombuffer(mp_failed_evts)
        NOE=np.frombuffer(mp_NOE)
        SX=np.frombuffer(mp_SX)
        SXI=np.frombuffer(mp_SXI)
        SXX=np.reshape(np.frombuffer(mp_SXX),(Speclength,Speclength))
        SI=np.frombuffer(mp_SI)
        SI2=np.frombuffer(mp_SI2)
        
        CovXX=((NOE-1)**-1)*(SXX-np.outer(SX,SX)*(NOE**-2))
   	
        Cov_Sliced += CovXX
    
        CovIX=(SXI-(SI*SX)/(NOE**2))/(NOE-1)
        CovII=(SI2-(SI*SI)/(NOE**2))/(NOE-1)
        pCov_Sliced += ((NOE-1)/(NOE-2))*(CovXX-1.05*np.outer(CovIX,CovIX)/CovII)

        NOE_tot += NOE
        SX_tot += SX
        SXI_tot += SXI
        SXX_tot += SXX
        SI_tot += SI
        SI2_tot += SI2

	#reset the values to 
        NOE = np.zeros(1,dtype=np.float64)
        SX = np.zeros(Speclength,dtype=np.float64)
        SXI = np.zeros(Speclength,dtype=np.float64)
        SXX = np.zeros((Speclength*Speclength),dtype=np.float64)
        SI = np.zeros(1,dtype=np.float64)
        SI2 = np.zeros(1,dtype=np.float64)
 
        print  'finished %i/%i files'%(nf+slicestep,len(fileslist))
    #pool.close()
    print 'failed events = %0.1f'%failed_evts, 'success events = %0.1f'%NOE_tot
    
    CovElements={'SX':SX_tot,'SXI':SXI_tot,'SXX':SXX_tot,'SI':SI_tot,'SI2':SI2_tot,'NOE':NOE_tot,'Cov':Cov_Sliced,'pCov':pCov_Sliced}
    
    np.savez(spath+'CovElements_%s_%ishotsx1000_Sliced%i.npz'%(fscan,nfilesstop,slicestep),**CovElements)
    print 'merge done'
    
    
    stop=time.time()-start
    print 'Completed in %0.3f' %stop 
    
#    try:    
#        processes=[mp.Process(target=CalculateCov,args=(argum,i)) for i in range(argum.cores)]
#        print processes
#    except:
#        print 'ERROR def processes'
#
#
#    for p in processes:
#        try:
#            p.start()
#            print p
#        except:
#            print 'ERROR launch processes'
#    for p in processes:
#        print p.is_alive()
#    for p in processes:
#        p.join()
#        
#    ## Merge the data all together    
#    try:
#        CalculateCovMerge(argum)
#    except:
#        print 'ERROR merging files'
#







