import numpy as np
import pandas as pd
import scipy.io as sio
import glob
import ctypes as c
import time


def myWorker(fname):
    filename=fname
    
    dframe=pd.read_table(filename)
    header=list(dframe.columns.values)[0]
    posX=[]
    posY=[]
    posT=[]

    counts=0
    eiCorr=[]
    eeCorr=[]
    iiCorr=[]
    
    
    print filename
    for S in dframe[header]:
    
        myshotarr=np.array(S.split(',')).astype(float)
        counts+=myshotarr[2]
        
        if myshotarr[2]>3:
            continue
        
        if myshotarr[2]>=2:
            T=myshotarr[4::4]
            mask=np.array((T>0.0)*(T<100.0))
            maskions=np.array((T>3200)*(T<15000))
            
            # Build electron-electron ion-ion correlations
            if len(T[maskions])==0:
                continue
            Xe=myshotarr[5::4][mask]
            Ye=myshotarr[6::4][mask]
            Re=np.sqrt(Xe**2+Ye**2)
            
            if len(Re)!=0:
                eeCorr.append([[i,j] for i in Re for j in Re])
            if len(T[~mask]!=0):
                iiCorr.append([[i,j] for i in T[~mask] for j in T[~mask]])
            
            if mask[0]==False:
                continue
            
            Xe=myshotarr[5]
            Ye=myshotarr[6]
            
            Re=np.sqrt(Xe**2+Ye**2)
            
            # Correlation between electron radius and ion TOF
            eiCorr.append([[Re,j] for j in T[1:]])
            
            posX.append([Xe]) # save electrons positions and ion TOF
            posY.append([Ye])
            posT.append(T[1:])
            
        elif myshotarr[2]==1:
            T=myshotarr[4]
            if (T>0.0)*(T<100.0)==False:
                continue
            
            X=myshotarr[5]
            Y=myshotarr[6]
            posX.append([X])
            posY.append([Y])
            posT.append([T])
            
        else:
            continue
        
    # Build the histograms
    
    HistVMIe,hvex,hvey=np.histogram2d(np.concatenate(posX),np.concatenate(posY),bins=(np.arange(-4000,4008,8),np.arange(-4000,4008,8)))
    hvcx=(hvex[1:]+hvex[:-1])/2.0
    hvcy=(hvey[1:]+hvey[:-1])/2.0
    
    HistPEPICO,hcex,hcey=np.histogram2d(np.concatenate(eiCorr).T[0],np.concatenate(eiCorr).T[1],bins=(np.arange(0,4008,8),np.arange(0,15004,4)))
    hccx=(hcex[1:]+hcex[:-1])/2.0
    hccy=(hcey[1:]+hcey[:-1])/2.0
    
    HistPEPECO,hcex,hcey=np.histogram2d(np.concatenate(eeCorr).T[0],np.concatenate(eeCorr).T[1],bins=(np.arange(0,4008,8),np.arange(0,4008,8)))
    hccx=(hcex[1:]+hcex[:-1])/2.0
    hccy=(hcey[1:]+hcey[:-1])/2.0
    
    HistPIPICO,hcex,hcey=np.histogram2d(np.concatenate(iiCorr).T[0],np.concatenate(iiCorr).T[1],bins=(np.arange(0,15004,4),np.arange(0,15004,4)))
    hccx=(hcex[1:]+hcex[:-1])/2.0
    hccy=(hcey[1:]+hcey[:-1])/2.0
    
    try:
        HistTOFi,htex=np.histogram(np.concatenate(posT),bins=(np.arange(0,15002,4)))
        htcx=(htex[1:]+htex[:-1])/2.0
    except:
        HistTOFi=np.zeros((3750))


    # return the histograms to the main function 
    
    return HistVMIe, HistPEPICO, HistPEPECO, HistPIPICO, HistTOFi


###....MAIN function of the program....###
if __name__ == '__main__':
    
    
    start=time.time()
    
    VMIe=np.zeros((1000,1000))
    PEPICO=np.zeros((500,3750))
    PEPECO=np.zeros((500,500))
    PIPICO=np.zeros((3750,3750))
    TOFi=np.zeros((3750))
    
    flist=glob.glob('/media/thomasb/New Volume/Data/EPFL/PEPICO/DLD/20180720_0830*.dat')[:]
    
    for f in flist:
        tmpVMIe,tmpPEPICO,tmpPEPECO,tmpPIPICO,tmpTOFi=myWorker(f)

        VMIe=tmpVMIe
        PEPICO+=tmpPEPICO
        PEPECO+=tmpPEPECO
        PIPICO+=tmpPIPICO
        TOFi+=tmpTOFi
        
    # Save all the maps in one file
    np.savez('/media/thomasb/New Volume/Data/EPFL/PEPICO/DLD/20180720_0830.npz',VMIe=VMIe,PEPICO=PEPICO,PEPECO=PEPECO,PIPICO=PIPICO,TOFi=TOFi)
    
    stop=time.time()-start
    print 'Completed in %0.3f' %stop 
    #print VMIe.shape
    
    

    
    
    
    
    
    
