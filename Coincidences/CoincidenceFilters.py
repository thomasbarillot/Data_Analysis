# -*- coding: utf-8 -*-
"""
Created on Fri Nov  4 10:13:33 2016

@author: tbarillot
"""

import numpy as np
import scipy.io as sio
import glob
import sys
import os

sys.path.append('/Users/tbarillot/Documents/Packages/igor/')
sys.path.append('/Users/tbarillot/Documents/PostDoc_imperial_LabMac/007_repositories/Coincidences/')

## Filter an array with other arrays conditions
def eFilter1D(arr,rawarr,flt_list,params_list):
    
    if len(flt_list)==0:
        print 'no filter specified'
        return None
    elif len(flt_list)==1:
        var=rawarr['%s'%flt_list[0]]
        params=params_list[0]
        arrflt=[arr[x] for x in range(0,len(arr)) if var[x]==params]
       
        return arrflt
    
    elif len(flt_list)==2:
        
        var1=rawarr['%s'%flt_list[0]]
        params1=params_list[0]
        var2=rawarr['%s'%flt_list[1]]
        params2=params_list[1]
        
        if len(params1)==1 and len(params2)==1:
            arrflt=[arr[x] for x in range(0,len(arr)) if var1[x]==params1 and var2[x,0]==params2]
        elif len(params1)==1 and len(params2)==2:
            arrflt=[arr[x] for x in range(0,len(arr)) if var1[x]==params1 and var2[x,0]>params2[0] and var2[x,0]<params2[1]]
        elif len(params1)==2 and len(params2)==2:
            arrflt=[arr[x] for x in range(0,len(arr)) if var1[x]>params1[0] and var1[x,0]<params1[1] and var2[x,0]>params2[0] and var2[x,0]<params2[1]]
        
        return arrflt



## Correct the radius based on a reference ratio   
def RadiusCorrection(arrR,arrPhi,corrf):
    
    idxarr=np.digitize(arrPhi,np.arange(-np.pi,np.pi,2*np.pi/180.0))
    arrRcorr=[arrR[x]/corrf[idxarr[x]-2] for x in range(0,len(arrR))]
    return arrRcorr
    
    
    
## calculate coincidence positions on coincidence map, autocorrelation value is either True or False
def MakeCoincidenceMap(arr1,arr2,autoco):
    MapX=[]
    MapY=[]
    for i in range(len(arr1[:,0])):
        subarr1=arr1[i,:]
        subarr2=arr2[i,:]
        if autoco:
            submap=[[[vk,vj] for k,vk in enumerate(subarr2) if k>=j] for j,vj in enumerate(subarr1)]
            submap=filter(None,submap)
        else:
            submap=[[[vk,vj] for k,vk in enumerate(subarr2) if k>j] for j,vj in enumerate(subarr1)]
            submap=filter(None,submap)
        MapX.append(np.concatenate(submap)[:,0])
        MapY.append(np.concatenate(submap)[:,1])
        
    MX=np.concatenate(MapX)
    MY=np.concatenate(MapY)
    
    return MX,MY

