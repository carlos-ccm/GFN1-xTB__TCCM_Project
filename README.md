# Semiempirical GFN1-xTB Method code

Development of a coding project for the TCCM Master. The project consist of creating from scratch a code that performs semiempirical DFT calculations for a given molecule. 
Author: Carlos Cruz Marin 
Mail: carloscruzmarin.bdn@gmail.com
X/Twitter: https://twitter.com/c_cruz_TC



## Table of Contents

1. [Installation](#installation)
2. [Inputs](#inputs)
3. [Usage](#usage)

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



