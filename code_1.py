#!/usr/bin/env python3

import sys 
import os
import numpy as np
import argparse
import time

def convert_fortran_to_python(num_string):
    return float(num_string.replace('d','e'))

def parameters_assigment(all_params):
    
    #enumerate gives us the chance of getting the index of each element in the list
    for i, line in enumerate(all_params): 
        if 'kf' in line:
            print(line)
            
            #When there is a match, we look for the following line for storing the values
            kf = float(all_params[i+1])
            print(kf)
            
        elif 'alpha' in line:
            print(line)
            alpha = all_params[i+1].split()
            
        elif 'Zeff' in line:
            print(line)
            zeff = all_params[i+1].split()
            
    return(kf,alpha,zeff)


def repulsion_energy(num_atoms):
    
    #Initialization of the repulsion energy term
    repulsion_energy = 0
    for a in range(num_atoms):
        for b in range(a+1, num_atoms):
            dist_ab = 0
            for i in range(3):
                dist_ab += (coordinates[a,i]-coordinates[b,i])**2
                
            dist_ab = np.sqrt(dist_ab)
            print("Distance between",atom_list[a],"and",atom_list[b],"is:",dist_ab)    
            #repulsion_energy += (z_dict.get(element_list[a]) * z_dict.get(element_list[b])/dist_ab) * np.exp(-np.sqrt() * dist_ab ** kf) 
    
    return(repulsion_energy)

if __name__ == '__main__': 
    st = time.time()
    
    #Input control 
    parser = argparse.ArgumentParser(description = 'xTB semiempirical DFT method for small molecules')
    parser.add_argument('-i', required= True, help='molecule coordinates input')
    parser.add_argument('-p', required=True, help='parameters of the elements')

    args = parser.parse_args()
    molecule = args.i
    parameters = args.p
    atom_list = []

    # Files reading
    with open(parameters,'r') as file:
       all_params = file.read().split('\n')
       print(all_params)
       
    with open(molecule,'r') as file:
        num_atoms = int(file.readline().strip())
        molecule_name = file.readline().strip()
        coordinates = np.zeros((num_atoms,3))
        
        for n in range(num_atoms):
           data = file.readline().strip()
           atom_list.append(data.split()[0])
           
           # After the first two lines, the atoms
           for i in range(3):
               coordinates[n][i] = data.split()[i+1]
           
    print("Number of atoms in the system:",num_atoms)
    print('Molecule of study:',molecule_name)
    for i in range(num_atoms):
        print("Atom of the system:",atom_list[i],"with coordinates:",coordinates[i,:])
    
    parameters_assigment(all_params)
    repulsion_energy(num_atoms)
    
    
    et = time.time()
    
    print("Execution time in seconds {:.2f}".format(et-st))
        
    

    
        
    

    