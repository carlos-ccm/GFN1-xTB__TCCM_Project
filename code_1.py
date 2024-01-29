#!/usr/bin/env python3

import numpy as np
import argparse
import time

##################################### TOOLS #####################################
def convert_fortran_to_python(num_string):
    return float(num_string.replace('d','e'))
###############################################################################################################

##################################### PARAMETERS FROM FILE #####################################
def parameters_assigment(all_params):
    alpha_dict,zeff_dict,q_dict,r_dict,ref_coord_num_dict,na_nb_dict,kll_dict,kcn_dict,hamiltonian_parameters_dict,electronegativities_dict,r_dict_2 = {},{},{},{},{},{},{},{},{},{},{}
    
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
        elif 'from eV to hartree' in line:
            energy_conversion = float(convert_fortran_to_python(all_params[i+1]))
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
        
        elif 'Unscaled covalent radii (in Angstrom) - for eq. 14 - ' in line:
            radii = all_params[i+1].split()
            for key, value in zip(element_list,radii):
                r_dict [key] = convert_fortran_to_python(value) * (4/3) * (1/distance_conversion)
    
        elif 'Reference coordination' in line:
            j = 0
            for key in element_list:
                j += 1   
                ref_coord_num_dict [key] = all_params[i+j].split() #This dict contains the reference coordination numbers CN_i^A per atom type
        
        elif 'Number of reference coordination' in line:
            ref_coord_num = all_params[i+1].split()
            for key,value in zip(element_list,ref_coord_num):
                na_nb_dict [key] = value #This dict contains the number of reference coordination numbers
        elif 'KAB term' in line:
            khh = (convert_fortran_to_python((all_params[i+1]).split()[0]))
            khn = (convert_fortran_to_python((all_params[i+1]).split()[1]))
        elif "kll'" in line:
            off_diagonal_kll_elements = int(all_params[i+1].split()[0])

            for j in range(2,off_diagonal_kll_elements+2):
                key1 = str((all_params[i+j]).split()[0])
                key2 = str((all_params[i+j]).split()[1])
                #initialize of the basis functions dictionary
                if key1 not in kll_dict:
                    kll_dict[key1] = {}

                #Ensure the inner dictionary has the basis_function_type key
                if key2 not in kll_dict[key1]:
                    kll_dict[key1][key2] = {} 
                kll_dict [key1][key2] = convert_fortran_to_python((all_params[i+j]).split()[2])
        elif 'Scaling constants' in line:
            kcn_dict['1'] = convert_fortran_to_python((all_params[i+1].split()[0]))
            kcn_dict['2'] = convert_fortran_to_python((all_params[i+1].split()[1]))
            kcn_dict['3'] = convert_fortran_to_python((all_params[i+1].split()[2]))
        elif 'values of eta_A and kappa' in line:
            for j in range(1,len(dictionary_atom_types)+1):
                key_atom_type = str((all_params[i+j]).split()[0])
                key_shell_type = str((all_params[i+j]).split()[1])

                 #initialize of the basis functions dictionary
                if key_atom_type not in hamiltonian_parameters_dict:
                    hamiltonian_parameters_dict[key_atom_type] = {}

                #Ensure the inner dictionary has the basis_function_type key
                if key2 not in hamiltonian_parameters_dict[key_atom_type]:
                    hamiltonian_parameters_dict[key_atom_type][key_shell_type] = {} 
                    
                hamiltonian_parameters_dict[key_atom_type][key_shell_type] = {
                    'n0': int(all_params[i+j].split()[2]),
                    'eta_Al': float(all_params[i+j].split()[3]),
                    'HAl': convert_fortran_to_python((all_params[i+j]).split()[4]) * (1/energy_conversion),
                    'kpoly' : convert_fortran_to_python((all_params[i+j]).split()[5])}      
        elif 'kEN' in line:
            kEN = convert_fortran_to_python(all_params[i+1])
        elif 'Electronegativities' in line:
            electronegativities = all_params[i+1].split()
            for key,value in zip(element_list,electronegativities):
                electronegativities_dict [key] = convert_fortran_to_python(value)   
        elif 'Unscaled covalent radii (in Angstrom) R^\Pi_A,cov. For eq. 27.' in line:
            radii_2 = all_params[i+1].split()
            for key, value in zip(element_list,radii_2):
                r_dict_2 [key] = convert_fortran_to_python(value) * (1/distance_conversion)
                     
    return(distance_conversion,energy_conversion,kf,alpha_dict,zeff_dict,q_dict,a1,a2,s6,s8,kcn,kl,r_dict,ref_coord_num_dict,na_nb_dict,khh,khn,kll_dict,kcn_dict,hamiltonian_parameters_dict,kEN,electronegativities_dict,r_dict_2)
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
            repulsion_energy += (convert_fortran_to_python(zeff_dict.get(atom_list[a])) * convert_fortran_to_python( zeff_dict.get(atom_list[b])))/(dist_ab) * np.exp(-np.sqrt(convert_fortran_to_python(alpha_dict.get(atom_list[a]))*convert_fortran_to_python(alpha_dict.get(atom_list[b]))) * dist_ab ** kf) 
    
    print("\n----------------------------------------------------\n| Repulsion energy is:",repulsion_energy,'Hartree |\n----------------------------------------------------\n')
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
        l_ij_acc = 0
        for atom_a in range(int(na_nb_dict.get(atom_list[a]))):
            for atom_b in range(int(na_nb_dict.get(atom_list[b]))):
                l_ij = np.exp(-kl * ((coord_number - convert_fortran_to_python(ref_coord_num_dict.get(atom_list[a])[atom_a]))**2)) 
                value = 1
                l_ij_acc += l_ij
       
    elif n == 8:
        c6 =   1     
        value = 3 * c6 * np.sqrt(convert_fortran_to_python(q_dict.get(atom_list[a])) * convert_fortran_to_python(q_dict.get(atom_list[b])))
    
    #This damping constant is just the rAB,0
    damping_constant = (9 * convert_fortran_to_python(q_dict.get(atom_list[a])) * convert_fortran_to_python(q_dict.get(atom_list[b]))) ** (1/4)
    
    return(value,damping_constant)

def damping_func(dist_ab,n,damping_constant):
    damp = (dist_ab**n)/((dist_ab**n) + (a1*damping_constant + a2) ** n)
    return(damp)

def coord_calc(a,num_atoms):
    coord_number = 0
    for b in range(num_atoms):
        dist_ab = 0
        if b != a:
            for i in range(3):
                dist_ab += (coordinates_bohr[a,i] - coordinates_bohr[b,i])**2
                
            dist_ab = np.sqrt(dist_ab)
            coord_number += (1 + np.exp(-kcn * ((r_dict.get(atom_list[a]) + r_dict.get(atom_list[b])) / dist_ab - 1)))** -1 
    return(coord_number)

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

def basis_set_select(all_params):
    
    match = False
    basis_functions_dict = {}

    
    #This part is for extracting specifically the part of parameters.dat that contains the information about the basis set
    basis_sets_info = []
    for line in all_params:
        if '(c) d contraction' in line:
            match = True
            basis_sets_info.append(line)
            
        elif 'Hamiltonian' in line:
            match = False
            break
        
        elif match is True:
            basis_sets_info.append(line)
    
    #Now we iterate again over the extracted part to obtain the basis set for each atom
    for i in range(1,number_of_atom_types * 2*3,3): #the number of basis functions is just the number of atom types x basis functions per atom type x number of lines describing each basis function  
        atom_type_number = basis_sets_info[i].split()[0]
        
        atom_type = basis_sets_info[i].split()[0]
        basis_function_type = basis_sets_info[i].split()[1]
        number_of_primitives = basis_sets_info[i].split()[2]
        
        basis_function_zeta = basis_sets_info[i+1].split()
        basis_function_contraction = basis_sets_info[i+2].split()
        
        #initialize of the basis functions dictionary
        if atom_type not in basis_functions_dict:
            basis_functions_dict[atom_type] = {}

        #Ensure the inner dictionary has the basis_function_type key
        if basis_function_type not in basis_functions_dict[atom_type]:
            basis_functions_dict[atom_type][basis_function_type] = {}

        #Update the inner dictionary with the extracted information
        basis_functions_dict[atom_type][basis_function_type] = {
            'Contraction': basis_function_contraction,
            'Exponent': basis_function_zeta}

    basis_functions = {}
    n = 0
    for i,at in enumerate(atom_list):
        for j in range(basis_functions_per_atom_type.get(at)):
            
            if j == 0:
                function_type = 0
            else:
                function_type = 1
            shell = list(basis_functions_dict[str(dictionary_atom_types.get(at))].items())[function_type][0]

            basis_functions[n] = {
            'Atom':at,
            'Atom index':i,
            'Function type':shell,
            'Contraction':basis_functions_dict[str(dictionary_atom_types.get(at))][str(shell)]['Contraction'],
            'Exponent':basis_functions_dict[str(dictionary_atom_types.get(at))][str(shell)]['Exponent']            
            } 
            n += 1

    return (basis_functions_dict,basis_functions)

def overlap_eval():
    overlap_matrix = np.zeros((num_basis_funcs, num_basis_funcs))
    key_pp, key_sp, key_ss, key_ps, key_p,key_s = False,False,False,False,False,False
    μ = 0
    v = 0
    while μ < num_basis_funcs:
        while v < num_basis_funcs:          
            dist_ab = 0
            a = basis_functions[μ]['Atom index']
            b = basis_functions[v]['Atom index']
            
            if basis_functions[μ]['Function type'][0] == '3' and basis_functions[v]['Function type'][0] == '3':
                key_pp = True
                key_p = True

            elif basis_functions[μ]['Function type'][0] == '3' and basis_functions[v]['Function type'][0] != '3':
                key_ps = True

            elif basis_functions[μ]['Function type'][0] != '3' and basis_functions[v]['Function type'][0] == '3':
                key_sp = True
                key_s = True

            elif basis_functions[μ]['Function type'][0] != '3' and basis_functions[v]['Function type'][0] != '3':
                key_ss = True
            
            for i in range(3):
                dist_ab += (coordinates_bohr[a,i]-coordinates_bohr[b,i])**2                            
            dist_ab = np.sqrt(dist_ab)
            
            for k in range(len(basis_functions[μ]['Contraction'])):
                for l in range(len(basis_functions[v]['Contraction'])):
                    #k and l iterate over all the primitives of each basis function μ and ν
                    #We need these terms for any integral, so we compute them now                     
                    ζ_k = float(basis_functions[μ]['Exponent'][k])
                    ζ_l = float(basis_functions[v]['Exponent'][l])
                    d_k = float(basis_functions[μ]['Contraction'][k])
                    d_l = float(basis_functions[v]['Contraction'][l])
                    
                    ζ = ζ_k + ζ_l 
                    ξ = (ζ_k * ζ_l)/ζ
                    
                    r_p = (ζ_k * coordinates_bohr[a,:] + ζ_l * coordinates_bohr[b,:]) / ζ
                    S_kμlv = np.exp(-ξ * (abs(dist_ab) ** 2)) * (np.pi/ζ) ** (3/2)   
                                  
                    if key_sp:
                       overlap_matrix[μ,v] += d_k * d_l *  (r_p[0] - coordinates_bohr[b,0]) * S_kμlv
                       overlap_matrix[μ,v+1] += d_k * d_l *  (r_p[1] - coordinates_bohr[b,1]) * S_kμlv
                       overlap_matrix[μ,v+2] += d_k * d_l *  (r_p[2] - coordinates_bohr[b,2]) * S_kμlv
                       
                    elif key_pp:
                        overlap_matrix[μ,v] += d_k * d_l *  ((r_p[0] - coordinates_bohr[a,0])* (r_p[0] - coordinates_bohr[b,0]) * S_kμlv + 1/(2*ζ)*S_kμlv) #x|x
                        overlap_matrix[μ+1,v+1] += d_k * d_l *  ((r_p[1] - coordinates_bohr[a,1])* (r_p[1] - coordinates_bohr[b,1]) * S_kμlv + 1/(2*ζ)*S_kμlv) #y|y
                        overlap_matrix[μ+2,v+2] += d_k * d_l *  ((r_p[2] - coordinates_bohr[a,2])* (r_p[2] - coordinates_bohr[b,2]) * S_kμlv + 1/(2*ζ)*S_kμlv) #z|z
                                                
                        overlap_matrix[μ,v+1] += d_k * d_l *  (r_p[0] - coordinates_bohr[b,0]) * (r_p[1] - coordinates_bohr[b,1])* S_kμlv #x|y
                        overlap_matrix[μ+1,v] = overlap_matrix[μ,v+1] #y|x
                        overlap_matrix[μ,v+2] += d_k * d_l *  (r_p[0] - coordinates_bohr[b,0]) * (r_p[2] - coordinates_bohr[b,2])* S_kμlv #x|z
                        overlap_matrix[μ+2,v] = overlap_matrix[μ,v+2] #z|x
                        overlap_matrix[μ+1,v+2] += d_k * d_l *  (r_p[0] - coordinates_bohr[b,0]) * (r_p[1] - coordinates_bohr[b,1])* S_kμlv #y|z
                        overlap_matrix[μ+2,v+1] = overlap_matrix[μ,v+1] #z|y

                    elif key_ps:
                        overlap_matrix[μ,v] += d_k * d_l *  (r_p[0] + coordinates_bohr[a,0]) * S_kμlv
                        overlap_matrix[μ+1,v] += d_k * d_l *  (r_p[1] + coordinates_bohr[a,1]) * S_kμlv
                        overlap_matrix[μ+2,v] += d_k * d_l *  (r_p[2] + coordinates_bohr[a,2]) * S_kμlv
                        
                    elif key_ss:
                        overlap_matrix[μ,v] += d_k * d_l * S_kμlv
                        
            if key_pp:
                key_pp = False
                v = v + 3 
            elif key_ps:
                key_ps = False
                v = v + 1
            elif key_sp:
                key_sp = False
                v = v + 3
            elif key_ss:
                key_ss = False
                v = v + 1
        if key_p:
            μ = μ + 3
            key_p = False
        elif key_s:
            μ += 1
            key_s = False
        v = 0
                
    print(overlap_matrix)
    return(overlap_matrix)

def inverse_square_root_overlap(overlap_matrix):
    eigenvalues, eigenvectors = np.linalg.eig(overlap_matrix)
    L = eigenvectors
    Λ = np.diag(eigenvalues) #this is the diagonal matrix Λ = L^-1 * S * L (L^-1 is denoted in the code as L_inv) 
    L_inv = np.linalg.inv(L) # as this is unitary, transpose of L is equivalent to inverse of L
    Λ_inv_sqrt = np.diag(1 / np.sqrt(eigenvalues))
    S_inv_sqrt = L @ Λ_inv_sqrt @ L_inv

    return(S_inv_sqrt)

def scaling_function(μ,v,basis_functions):
    if  (basis_functions[μ]['Atom'] == 'H' and basis_functions[v]['Atom'] == 'H'):
        scaling_parameter = khh
    elif (basis_functions[μ]['Atom'] == 'N' and basis_functions[v]['Atom'] == 'H') or (basis_functions[μ]['Atom'] == 'H' and basis_functions[v]['Atom'] == 'N'):
        scaling_parameter = khn
    else:
        scaling_parameter = 1 
           
    return (scaling_parameter)

def electronic_energy(basis_functions):
    ##Calculation of the electronic term is shown below ##    
        #First the mu|H0|nu term, from the 0th order Hamiltonian
    H0 = np.zeros((num_basis_funcs,num_basis_funcs))
    for μ in range(num_basis_funcs):
       for v in range(num_basis_funcs):
            h_l_a = hamiltonian_parameters_dict[str(dictionary_atom_types.get(basis_functions[μ]['Atom']))][basis_functions[μ]['Function type'][0]]['HAl'] * (1 + kcn_dict[basis_functions[μ]['Function type'][0]] * coord_calc(basis_functions[μ]['Atom index'], num_atoms)) #effective atomic energy level h_l_a 
            h_l_b = hamiltonian_parameters_dict[str(dictionary_atom_types.get(basis_functions[v]['Atom']))][basis_functions[v]['Function type'][0]]['HAl'] * (1 + kcn_dict[basis_functions[v]['Function type'][0]] * coord_calc(basis_functions[v]['Atom index'], num_atoms)) #effective atomic energy level h_l_b
            if μ == v:
                H0[μ,v] = hamiltonian_parameters_dict[str(dictionary_atom_types.get(basis_functions[μ]['Atom']))][basis_functions[μ]['Function type'][0]]['HAl'] * (1 + kcn_dict[basis_functions[μ]['Function type'][0]] * coord_calc(basis_functions[μ]['Atom index'], num_atoms)) #effective atomic energy level h_l_a 
            else:
                scaling_parameter = scaling_function(μ,v,basis_functions)
                dist_ab = 0
                a = basis_functions[μ]['Atom index']
                b = basis_functions[v]['Atom index']
                for i in range(3):
                    dist_ab += (coordinates_bohr[a,i]-coordinates_bohr[b,i])**2                            
                dist_ab = np.sqrt(dist_ab)
                R_AB = r_dict_2.get(basis_functions[μ]['Atom']) + r_dict_2.get(basis_functions[v]['Atom'])
                        
                pi_rab_ll = (1 + hamiltonian_parameters_dict[str(dictionary_atom_types.get(basis_functions[μ]['Atom']))][basis_functions[μ]['Function type'][0]]['kpoly'] * ((dist_ab/R_AB)**0.5)) * (1 + hamiltonian_parameters_dict[str(dictionary_atom_types.get(basis_functions[v]['Atom']))][basis_functions[v]['Function type'][0]]['kpoly'] * ((dist_ab/R_AB)**0.5))
                if basis_functions[μ]['Function type'][0] == '2' or basis_functions[v]['Function type'][0] == '2':
                    variation_electronegativity_term = 1
                else:
                    variation_electronegativity_term = 1 + kEN * (electronegativities_dict.get(basis_functions[μ]['Atom']) - electronegativities_dict.get(basis_functions[v]['Atom']))**2
                
                if float(basis_functions[μ]['Function type'][0]) <= float(basis_functions[v]['Function type'][0]):
                    H0[μ,v]  = scaling_parameter * kll_dict[str(basis_functions[μ]['Function type'][0])][str(basis_functions[v]['Function type'][0])] * 0.5 * (h_l_a + h_l_b) * overlap_matrix[μ,v] * variation_electronegativity_term * pi_rab_ll
                else:
                    H0[μ,v]  = scaling_parameter * kll_dict[str(basis_functions[v]['Function type'][0])][str(basis_functions[μ]['Function type'][0])] * 0.5 * (h_l_a + h_l_b) * overlap_matrix[μ,v] * variation_electronegativity_term * pi_rab_ll                
    print("\n0th Hamiltonian\n")
    print(H0)
    return(H0)    
        
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
           atom_list.append(data.split()[0]) #atom list is the list of atoms of the system
           
           #After the first two lines, the atom's coordinates in Angstroms
           for i in range(3):
               coordinates[n][i] = data.split()[i+1]
           
    print("Number of atoms in the system:",num_atoms)
    print('Molecule of study:',molecule_name)
    
    
    element_list = ['H','C','N','O'] #list of elements present in the system
    number_of_atom_types = len(element_list)

    dictionary_atom_types = {'H':1,'C':2,'N':3,'O':4}
    electrons_per_atom_type = {'H':1,'C':4,'N':5,'O':6}
    basis_functions_per_atom_type = {'H':2,'C':4,'N':4,'O':4} #3 per p type
    
    system_elec = 0
    num_basis_funcs = 0
    
    for i in range(num_atoms):
        system_elec += electrons_per_atom_type.get(atom_list[i])
        num_basis_funcs += basis_functions_per_atom_type.get(atom_list[i])
        
    print('Number of electrons and occupied orbitals:',system_elec,int(system_elec/2))
    print('Total number of shells and basis functions:',num_atoms*2,num_basis_funcs)
    
    distance_conversion,energy_conversion,kf,alpha_dict,zeff_dict,q_dict,a1,a2,s6,s8,kcn,kl,r_dict,ref_coord_num_dict,na_nb_dict,khh,khn,kll_dict,kcn_dict,hamiltonian_parameters_dict,kEN,electronegativities_dict,r_dict_2= parameters_assigment(all_params)
    coordinates_bohr = coordinates * 1 / distance_conversion
    
    for i in range(num_atoms):
        print("Atom of the system:",atom_list[i],"with coordinates:",coordinates_bohr[i,:],'in bohrs')
    
    #Simple terms of the energy, repulsion and dispersion
    repulsion_contribution()
    dispersion_contribution()
    
    
    #Now we calculate the components of the electronic energy
    basis_functions_dict, basis_functions = basis_set_select(all_params)
    overlap_matrix = overlap_eval()
    S_inv_sqrt = inverse_square_root_overlap(overlap_matrix)
    
    H0 = electronic_energy(basis_functions)
    
    
    et = time.time()
    
    print("Execution time in seconds {:.2f}".format(et-st))  