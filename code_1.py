#!/usr/bin/env python3

import sys 
import os
import numpy as np
import argparse
import time

if __name__ == '__main__': 
    st = time.time()
    
    #Input control (NOTE:We will not use mol2 files)
    parser = argparse.ArgumentParser(description = 'xTB semiempirical DFT method for small molecules')
    parser.add_argument('-i', required= True, help='molecule coordinates input')
    parser.add_argument('-p', required=True, help='parameters of the elements')

    args = parser.parse_args()
    test = args.i

    if test.endswith('.sdf') or test.endswith('.sd'):
        suppl = Chem.SDMolSupplier(test, removeHs=False)

    model_dir = os.getcwd()+'/'