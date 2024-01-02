#!/usr/bin/env python3

import sys 
import os
import numpy as np
import argparse
import time

def convert_fortran_to_python(num_string):
    return float(num_string.replace('d','e'))

##################################### PARAMETERS FROM FILE #####################################
def parameters_assigment(all_params):
    alpha_dict = {}
    zeff_dict = {}
    q_dict = {}
    r_dict = {}
    
    #enumerate gives us the chance of getting the index of each element in the list
    for i, line in enumerate(all_params): 
        if 'kf' in line:
            #When there is a match, we look for the following line for storing the values in dictionaries using the function zip
            kf = float(all_params[i+1])
            
        elif 'alpha' in line: 
            alpha = all_params[i+1].split()
            for key,value in zip(element_list,alpha):
                alpha_dict [key] = value
            
        elif 'Zeff' in line:  
            zeff = all_params[i+1].split()
            for key,value in zip(element_list,zeff):
                zeff_dict [key] = value
                
        elif 'bohr to Angstrom' in line:
            distance_conversion = float(convert_fortran_to_python(all_params[i+1]))
        
        elif 'QA values' in line:
            q = all_params[i+1].split()
            for key, value in zip(element_list,q):
                q_dict [key] = value
        
        elif 'Overall parameters' in line:
            a1 = (convert_fortran_to_python((all_params[i+1]).split()[0]))
            a2 = (convert_fortran_to_python((all_params[i+1]).split()[1]))
            s6 = (convert_fortran_to_python((all_params[i+1]).split()[2]))
            s8 = (convert_fortran_to_python((all_params[i+1]).split()[3]))
            kcn = (convert_fortran_to_python((all_params[i+1]).split()[4]))
            kl = (convert_fortran_to_python((all_params[i+1]).split()[5]))
        
        elif 'covalent radii' in line:
            radii = all_params[i+1].split()
            for key, value in zip(element_list,radii):
                r_dict [key] = convert_fortran_to_python(value) * (4/3) * (1/distance_conversion)
    

    return(distance_conversion,kf,alpha_dict,zeff_dict,q_dict,a1,a2,s6,s8,kcn,kl,r_dict)
###############################################################################################################

##################################### REPULSION TERM #####################################

def repulsion_contribution():
    
    #Initialization of the repulsion energy term
    repulsion_energy = 0
    for a in range(num_atoms):
        for b in range(a+1, num_atoms):
            dist_ab = 0
            for i in range(3):
                dist_ab += (coordinates_bohr[a,i]-coordinates_bohr[b,i])**2
                
            dist_ab = np.sqrt(dist_ab)
            print("Distance between",atom_list[a],"and",atom_list[b],"is:",dist_ab)    
            repulsion_energy += (convert_fortran_to_python(zeff_dict.get(atom_list[a])) * convert_fortran_to_python( zeff_dict.get(atom_list[b])))/(dist_ab) * np.exp(-np.sqrt(convert_fortran_to_python(alpha_dict.get(atom_list[a]))*convert_fortran_to_python(alpha_dict.get(atom_list[b]))) * dist_ab ** kf) 
    
    print("Repulsion energy is:",repulsion_energy,'Hartree')
    return(repulsion_energy)
###############################################################################################################

##################################### DISPERSION TERM #####################################

def switcher(n):
    if n == 6:
        value = s6
    elif n == 8:
        value = s8
    return(value)

def cn(a,b,n,coord_number):

    if n == 6:
        value = 1
    elif n == 8:
        c6 =   1     
        value = 3 * c6 * np.sqrt(convert_fortran_to_python(q_dict.get(atom_list[a])) * convert_fortran_to_python(q_dict.get(atom_list[b])))
    
    damping_constant = (9 * convert_fortran_to_python(q_dict.get(atom_list[a])) * convert_fortran_to_python(q_dict.get(atom_list[b]))) ** (1/4)
    
    return(value,damping_constant)

def damping_func(dist_ab,n,damping_constant):
    
    damp = (dist_ab**n)/((dist_ab**n) + (a1*damping_constant + a2) ** n)
    return(damp)

def coord_calc(a,num_atoms):
    coord_number = 0
    for b in range(num_atoms):
        if b != a:
            coord_number += (1 + np.exp(-kcn * r_dict.get(atom_list[a]) + r_dict.get(atom_list[b])) ** -1) 

def dispersion_contribution():
    
    #Initialization of the dispersion energy term
    dispersion_energy = 0
    for a in range(num_atoms):
        coord_number = coord_calc(a,num_atoms)
        for b in range (a+1, num_atoms):
            dist_ab = 0
            for i in range(3):
                dist_ab += (coordinates_bohr[a,i]-coordinates_bohr[b,i])**2
                
            dist_ab = np.sqrt(dist_ab)
            for n in [6, 8]:
                cn_n,damping_constant= cn(a,b,n,coord_number)
                dispersion_energy += ((-switcher(n) * cn_n)/(dist_ab**n)) * damping_func(dist_ab,n,damping_constant)
            
    
    return(dispersion_energy)
###############################################################################################################

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
       
          
       
    with open(molecule,'r') as file:
        num_atoms = int(file.readline().strip())
        molecule_name = file.readline().strip()
        coordinates = np.zeros((num_atoms,3))
        
        for n in range(num_atoms):
           data = file.readline().strip()
           atom_list.append(data.split()[0])
           
           # After the first two lines, the atom's coordinates in Angstroms
           for i in range(3):
               coordinates[n][i] = data.split()[i+1]
           
    print("Number of atoms in the system:",num_atoms)
    print('Molecule of study:',molecule_name)
    for i in range(num_atoms):
        print("Atom of the system:",atom_list[i],"with coordinates:",coordinates[i,:])
    
    element_list = ['H','C','N','O']
    
    distance_conversion,kf,alpha_dict,zeff_dict,q_dict,a1,a2,s6,s8,kcn,kl,r_dict= parameters_assigment(all_params)
    coordinates_bohr = coordinates * 1 / distance_conversion
    
    repulsion_contribution()
    dispersion_contribution()
    
    
    
    et = time.time()
    
    print("Execution time in seconds {:.2f}".format(et-st))  