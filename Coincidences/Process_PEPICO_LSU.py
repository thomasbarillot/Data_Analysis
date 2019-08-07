import numpy as np
import pandas as pd
import scipy.io as sio

import ctypes as c

## mp imports
import multiprocessing as mp
from multiprocessing import Pool

## arguments parser
import argparse 
parser = argparse.ArgumentParser()
parser.add_argument("--fpath", help="enter filepath",type=int)
args=parser.parse_args()






if __name__ == '__main__':
    
    pool=Pool(processes=8)
    dpath=args.fpath

    dframe=pd.read_table(dpath)
    header=list(dframe.columns.values)[0]
    
    VMIelectrons_arr=mp.Array(c.c_long,1000*1000)
    VMIions_arr=mp.Array(c.c_long,1000*1000)
    
    
    
    
    
    
