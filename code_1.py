#!/usr/bin/env python3

import sys 
import os
import numpy as np
import argparse
import time

if __name__ == '__main__': 
    st = time.time()
    
    #Input control (NOTE:We will not use mol2 files)
    parser = argparse.ArgumentParser(description = 'Prediction of atomic logP using atomic fingerprint')
    parser.add_argument('-i', required= True, help='train_input file')

    args = parser.parse_args()
    test = args.i

    if test.endswith('.sdf') or test.endswith('.sd'):
        suppl = Chem.SDMolSupplier(test, removeHs=False)

    model_dir = os.getcwd()+'/'