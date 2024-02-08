# Semiempirical GFN1-xTB Method code

Development of a coding project for the TCCM Master. The project consist of creating from scratch a code that performs semiempirical DFT calculations for a given molecule, composed by H, C, N and O. Method based on the theory developed in “Extended tight-binding quantum
chemistry methods”, C. Bannwarth, E. Caldeweyher, S. Ehlert, A. Hansen, P. Pracht, J. Seibert,S. Spicher and S. Grimme, WIREs Comput. Mol. Sci. 2021, 11 e1493, DOI: https://doi.org/10.1002/wcms.1493.

Author: Carlos Cruz Marin

Mail: carloscruzmarin.bdn@gmail.com

X/Twitter: https://twitter.com/c_cruz_TC



## Table of Contents

1. [Installation](#installation)
2. [Inputs](#inputs)
3. [Usage](#usage)
4. [References](#references)

## Installation

For running this code you will only need the following packages installed:
- Numpy


```bash

# Set up steps
git clone https://github.com/carlos-ccm/xtb.git
pip install numpy

```


## Inputs

When running the code you will need to provide some inputs.
- Molecule geometry (.xyz format)
- File containing the parameters (parameters.dat)

## Usage

The code includes the following flags
```python

-help Description of the code
-i <xyz file>
-p <parameters.dat>
```

##References

-  A general overview of the GFNn-xTB methods (n = 0, 1, 2): “Extended tight-binding quantum
chemistry methods”, C. Bannwarth, E. Caldeweyher, S. Ehlert, A. Hansen, P. Pracht, J. Seibert,
S. Spicher and S. Grimme, WIREs Comput. Mol. Sci. 2021, 11 e1493, DOI: https://doi.org/
10.1002/wcms.1493.
-  What I refer to above as the ‘original paper’, which contains a systematic description of the
GFN1-xTB method and its performance: “A Robust and Accurate Tight-Binding Quantum
Chemical Method for Structures, Vibrational Frequencies, and Noncovalent Interactions of Large
Molecular Systems Parametrized for All spd-Block Elements (Z = 1 – 86)”, S. Grimme, C.
Bannwarth and P. Shushkov, J. Chem. Theory Comput. 2017, 13, 1989 – 2009, https://doi.
org/10.1021/acs.jctc.7b00118.
-  Description of the DFTD3 dispersion correction: “A consistent and accurate ab initio parametrization of density functional dispersion correction (DFT-D) for the 94 elements H-Pu”, S. Grimme,
J. Antony, S. Ehrlich and H. Krieg, J. Chem. Phys., 2010, 132, 154104, https://doi.org/10.
1063/1.3382344.
-  Description of the ‘Becke-Johnson’ short-range dampiung function: “A post-Hartree-Fock model
of intermolecular interactions: Inclusion of higher-order corrections”, E. R. Johnson and A. D.
Becke, J. Chem. Phys., 2006, 124, 174104, https://doi.org/10.1063/1.2190220,andrefs.
\therein..
-   A method for evaluating integrals over basis functions (only the simplest case of basis function
overlap integrals is needed for this project): “E


