
from matplotlib.colors import LinearSegmentedColormap,Normalize
import numpy as np
from numpy import ma
import matplotlib.cbook as cbook

cdictalps1 = {'red': ((0.0, 1.0, 1.0),
        (0.16667, 0.0, 0.0),
        (0.33333, 0.5, 0.5),
        (0.5, 0.0, 0.0),
        (0.66667, 1.0, 1.0),
        (1, 1.0, 1.0)),
    'green': ((0.0, 0.0, 0.0),
        (0.16667, 0.0, 0.0),
        (0.33333, 1.0, 1.0),
        (0.5, 0.5, 0.5),
        (0.66667, 1.0, 1.0),
        (0.83333, 0.0, 0.0),
        (1.0, 1.0, 1.0)),
    'blue': ((0.0, 1.0, 1.0),
        (0.33333, 1.0, 1.0),
        (0.5, 0.0, 0.0),
        (0.83333, 0.0, 0.0),
        (1.0, 1.0, 1.0))}

Alps1 = LinearSegmentedColormap('Alps1', cdictalps1)

cdictalps1pos = {'red': ((0.0, 0.0, 0.0),
        (0.33333, 1.0, 1.0),
        (1.0, 1.0, 1.0)),
    'green': ((0.0, 0.5, 0.5),
        (0.33333, 1.0, 1.0),
        (0.66667, 0.0, 0.0),
        (1.0, 1.0, 1.0)),
    'blue': ((0.0, 0.0, 0.0),
        (0.66667, 0.0, 0.0),
        (1.0, 1.0, 1.0))}
        


Alps1Pos = LinearSegmentedColormap('Alps1Pos', cdictalps1pos)

cdictalps1neg = {'red': ((0.0, 1.0, 1.0),
        (0.33333, 0.0, 0.0),
        (0.66667, 0.5, 0.5),
        (1, 0.0, 0.0)),
    'green': ((0.0, 0.0, 0.0),
        (0.33333, 0.0, 0.0),
        (0.66667, 1.0, 1.0),
        (1, 0.5, 0.5)),
    'blue': ((0.0, 1.0, 1.0),
        (0.66667, 1.0, 1.0),
        (1, 0.0, 0.0))}

Alps1Neg = LinearSegmentedColormap('Alps1Neg', cdictalps1neg)

class CovarianceNorm(Normalize):

    def __init__(self, vmin=None, vmax=None, clip=False,fpos=(lambda x: x**0.5),finvpos=(lambda x: x**2),fneg=None,finvneg=None):
        if fneg is None:
            fneg=fpos
        if finvneg is None:
            finvneg=finvpos
        if vmin is not None and vmax is not None:
            if vmin > vmax:
                raise ValueError("vmin must be less than vmax")
        self.fneg=fneg
        self.fpos=fpos
        self.finvneg=finvneg
        self.finvpos=finvpos
        Normalize.__init__(self, vmin, vmax, clip)
               
    """
    Normalize a given value to the 0-1 range on a log scale
    """
    def __call__(self, value, clip=None):
        if clip is None:
            clip = self.clip
         
        result, is_scalar = self.process_value(value) 
        
        vmin = self.vmin 
        vmax = self.vmax
        
        
        result[result>vmax]=vmax
        result[result<vmin]=vmin
        
        maskpositive=result>0
        masknegative=result<0
        if vmax>0 and vmin<0:                   
            result[masknegative]=-self.fneg(result[masknegative]/vmin)      
            result[maskpositive]=self.fpos(result[maskpositive]/vmax)
            
        elif vmax>0 and vmin>=0:                   
            result[maskpositive]=self.fpos((result[maskpositive]-vmin)/(vmax-vmin))     

        elif vmax<=0 and vmin<0:                   
            result[masknegative]=-self.fneg((result[maskpositive]-vmax)/(vmin-vmax))    
        
        result=(result+1)/2        
        
        self.autoscale_None(result)        
        return result

    def inverse(self, value):

        
        value=value*2-1
        vmin = self.vmin 
        vmax = self.vmax
           
        if cbook.iterable(value):
            
            
            maskpositive=value>0
            masknegative=value<0
           
            if vmax>0 and vmin<0:                   
                value[masknegative]=self.finvneg(-value[masknegative])*vmin    
                value[maskpositive]=self.finvpos(value[maskpositive])*vmax
            
            elif vmax>0 and vmin>=0:                   
                value[maskpositive]=self.finvpos(value[maskpositive])*(vmax-vmin)+vmin  
                value[masknegative]=-self.finvpos(value[masknegative])*(vmax-vmin)+vmin                 
            elif vmax<=0 and vmin<0:                   
                value[masknegative]=self.finvneg(-value[masknegative])*(vmin-vmax)+vmax

        else:
    
            if vmax>0 and vmin<0:     
                if value<0:
                    value=self.finvneg(-value)*vmin   
                else:
                    value=self.finvpos(value)*vmax
            
            elif vmax>0 and vmin>=0:            
                if value>0:            
                    value=self.finvpos(value)*(vmax-vmin)+vmin    
                    
            elif vmax<=0 and vmin<0:        
                if value<0:    
                    value=self.finvneg(-value)*(vmin-vmax)+vmax               
        return value
        
    def ticks(self,N=11):
        return self.inverse(np.linspace(0,1,N))

    def autoscale(self, A): 
        """
        vmin = self.vmin 
        vmax = self.vmax
        
        if vmax==0 or ma.max(A)==0:
            self.vmin = ma.min(A)
            self.vmax = -self.vmin
        elif vmin==0 or ma.min(A)==0:          
            self.vmax = ma.max(A) 
            self.vmin = -self.vmax
        
        else:    
            self.vmin = ma.min(A)
            self.vmax = ma.max(A)  
        """
        self.vmin = ma.min(A)
        self.vmax = ma.max(A)                 
        
        
    def autoscale_None(self, A):
        if self.vmin is not None and self.vmax is not None:
            return
        if self.vmin is None:
            self.vmin = ma.min(A)
        if self.vmax is None:
            self.vmax = ma.max(A)
            
            
class PositiveCovarianceNorm(Normalize):

    def __init__(self, vmin=None, vmax=None, clip=False,fpos=(lambda x: x**0.5),finvpos=(lambda x: x**2)):
        if vmin is not None and vmax is not None:
            if vmin > vmax:
                raise ValueError("vmin must be less than vmax")
        self.fpos=fpos
        self.finvpos=finvpos
        Normalize.__init__(self, vmin, vmax, clip)
               
    def __call__(self, value, clip=None):
        if clip is None:
            clip = self.clip
         
        result, is_scalar = self.process_value(value) 
        
        vmin = self.vmin 
        vmax = self.vmax
               
        result[result>vmax]=vmax
        result[result<vmin]=vmin
                
        result=self.fpos((result-vmin)/(vmax-vmin)) 
                        
        self.autoscale_None(result)        
        return result

    def inverse(self, value):
       
        vmin = self.vmin 
        vmax = self.vmax
        
        if cbook.iterable(value):            
            value=self.finvpos(value)*(vmax-vmin)+vmin            
        else:
            value=self.finvpos(value)*(vmax-vmin)+vmin             
        return value
        
    def ticks(self,N=11):
        return self.inverse(np.linspace(0,1,N))

    def autoscale(self, A): 
        self.vmin = ma.min(A)
        self.vmax = ma.max(A)                 
        
        
    def autoscale_None(self, A):
        if self.vmin is not None and self.vmax is not None:
            return
        if self.vmin is None:
            self.vmin = ma.min(A)
        if self.vmax is None:
            self.vmax = ma.max(A)
            
            
class NegativeCovarianceNorm(Normalize):

    def __init__(self, vmin=None, vmax=None, clip=False,fneg=(lambda x: x**0.5),finvneg=(lambda x: x**2)):

        if vmin is not None and vmax is not None:
            if vmin > vmax:
                raise ValueError("vmin must be less than vmax")
        self.fneg=fneg
        self.finvneg=finvneg
        Normalize.__init__(self, vmin, vmax, clip)
               
    def __call__(self, value, clip=None):
        if clip is None:
            clip = self.clip
         
        result, is_scalar = self.process_value(value) 
        
        vmin = self.vmin 
        vmax = self.vmax
               
        result[result>vmax]=vmax
        result[result<vmin]=vmin
                       
        result=-self.fneg((result-vmax)/(vmin-vmax))         
        result=result+1
                        
        self.autoscale_None(result)        
        return result

    def inverse(self, value):
       
        vmin = self.vmin 
        vmax = self.vmax
        
        value=value-1
        
        
        if cbook.iterable(value):            
            value=self.finvneg(-value)*(vmin-vmax)+vmax       
        else:
            value=self.finvneg(value)*(vmin-vmax)+vmax  

        return value
        
    def ticks(self,N=11):
        return self.inverse(np.linspace(0,1,N))

    def autoscale(self, A): 
        self.vmin = ma.min(A)
        self.vmax = ma.max(A)                 
        
        
    def autoscale_None(self, A):
        if self.vmin is not None and self.vmax is not None:
            return
        if self.vmin is None:
            self.vmin = ma.min(A)
        if self.vmax is None:
            self.vmax = ma.max(A)
       
class CovarianceNormRoot(CovarianceNorm):

    def __init__(self, vmin=None, vmax=None, clip=False,orderpos=2,orderneg=None ):
        if orderneg is None:
            orderneg=orderpos
        CovarianceNorm.__init__(self, vmin, vmax, clip ,fneg=(lambda x: x**(1./orderneg)),finvneg=(lambda x: x**(orderneg)),fpos=(lambda x: x**(1./orderpos)),finvpos=(lambda x: x**(orderpos)))
        
class PositiveCovarianceNormRoot(PositiveCovarianceNorm):

    def __init__(self, vmin=None, vmax=None, clip=False,orderpos=2 ):

        PositiveCovarianceNorm.__init__(self, vmin, vmax, clip ,fpos=(lambda x: x**(1./orderpos)),finvpos=(lambda x: x**(orderpos)))
        
class NegativeCovarianceNormRoot(NegativeCovarianceNorm):

    def __init__(self, vmin=None, vmax=None, clip=False,orderneg=2 ):

        NegativeCovarianceNorm.__init__(self, vmin, vmax, clip ,fneg=(lambda x: x**(1./orderneg)),finvneg=(lambda x: x**(orderneg)))

        
