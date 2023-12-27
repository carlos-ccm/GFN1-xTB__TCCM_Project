#!/usr/bin/env python3

import sys 
import os
import numpy as np
import argparse
import time
import pandas


def repulsion_energy(num_atoms):
    for n in num_atoms:
        repulsion_energy = 0
    
    
    
    return(repulsion_energy)

if __name__ == '__main__': 
    st = time.time()
    
    #Input control (NOTE:We will not use mol2 files)
    parser = argparse.ArgumentParser(description = 'xTB semiempirical DFT method for small molecules')
    parser.add_argument('-i', required= True, help='molecule coordinates input')
    parser.add_argument('-p', required=True, help='parameters of the elements')

    args = parser.parse_args()
    molecule = args.i
    parameters = args.p

    # Files reading
    with open(parameters,'r') as file:
       all_params = file.read()
       
    with open(molecule,'r') as file:
        file.readline()
        molecule_name = file.readline().strip()

        
    # The number of atoms is always the first line of a xyz file
    num_atoms = int(np.loadtxt(molecule, max_rows = 1))
    
    # After the first 2 rows, then the coordinates of the system appear
    coordinates = np.loadtxt(molecule,skiprows=2,usecols=(1,3))
    
    
    print(num_atoms)
    print('Molecule of study:',molecule_name)
    

    