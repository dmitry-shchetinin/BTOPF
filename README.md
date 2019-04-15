# BTOPF: package for tightening bounds on variables for AC OPF problem

## Overview 
BTOPF is a Matlab-based package for efficiently producing tighter bounds on variables in the AC OPF problem. Two types of bounds 
are currently considered: bounds on voltage phase angle differences and bounds on voltage magnitude differences for all branches. 
The user can choose the type(s) of bounds to tighten and the method(s) to be used. This repository includes the following files:

1. Matlab scripts required to run the bounds tightening algorithm.
2. Source code for mex functions, which can significantly speed up parts of code.
3. Pre-built mex files for Windows (in 'release', coming soon).
4. Example scripts.
5. Installation script. It is required only for compiling mex functions from source code.

- - - -

## Basic Requirements
Since the basic version of the package relies only on standard Matlab functions, no installation is necessary. Note that BTOPF uses one 
function from [MATPOWER](http://www.pserc.cornell.edu/matpower/), namely `makeYbus`. Therefore, you should either have MATPOWER installed or put this function in some 
folder in the Matlab path. To be able to run BTOPF every time you start Matlab, add folder `BTOPF\matlab` to Matlab path and save it.


## Making BTOPF faster
Currently, there are two optional things you could do to accelerate BTOPF:

1. Use mex functions. Parts of the code that are slow in Matlab were converted to so-called mex functions. The source code for these 
functions is available from this repository. To compile them, run script `install_BTOPF` and follow prompts. This requires that Matlab 
have a C compiler. The list of supported compilers is available [here]( https://ch.mathworks.com/support/compilers.html). Alternatively,
you can download pre-built mex functions from this repository (they are in 'release') and put them in folder `BTOPF\matlab`.
2. Use KLU function, which is part of [SuiteSparse](http://faculty.cse.tamu.edu/davis/suitesparse.html) distribution. This function does 
the LU matrix factorization and is particularly efficient for matrices arising in circuit-related problems.

Using both of these options helps make BTOPF 2-4 times faster.

- - - -

## Bound tightening algorithm
The following bound tightening (BT) methods are available:

1. BT based on the feasible set of the thermal limit constraint. This is the fastest method and its complexity is linear with respect to 
the system size. It is only applicable to branches with known values of the thermal limit. The bounds on both angle and voltage differences 
can be tightened.
2. BT based on convex envelopes of power flows. The bounds on both angle and voltage differences can be tightened. Each bound is tightened by 
solving a linear optimization problem (LP), the size of which is ~6-8 times the number of buses in the system. The LP is formulated such that 
it admits a fast solution. The bounds on angle differences can be tightened iteratively. This is the slowest method.
3. BT based on convex envelopes of bus current injections. Only the bounds on voltage differences can be tightened. Each bound is tightened by 
solving an LP and/or second-order cone problem, the size of which is ~2 times the number of buses in the system. The problems are formulated such that 
they admit a fast solution.

These methods complement each other and the best result can be achieved when using all three. The description of algorithm's inputs and outputs 
is given below.

### Inputs
The inputs include the data of the test system and algorithm's options.

#### Test system (required)
Structure containing the system model in MATPOWER format.

#### Algorithm's options (optional)
Structure whose fields are algorithm's options. If this argument is not provided or left empty, the default options are used. The user only needs 
to pass the non-default options, the rest are initizialed to their default values. The following options are available:

- type of bounds to be tightened:
    - 1 - angle differences (default),
    - 2 - voltage differences,
    - 3 - both angle and voltage differences.
- BT method to be used:
    - 1 - all applicable methods (default),
    - 2 - BT based on line flow constraints,
    - 3 - BT based on convex envelopes of power flows,
    - 4 - BT based on envelopes of bus injections (only applicable for tightening angle bounds).
 - linear solver to be used:
    - 1 - default Matlab LU (default),
    - 2 - KLU from SuiteSparse (must be installed separately).
 - usage of mex functions to speed up computations:
    - 0 - mex functions are not used (default),
    - 1 - mex functions are used (must compile them locally or download pre-built binaries and put them into path).
 - number of iterations for tightening methods that are based on convex envelopes (default is 5).
 - stopping condition for iterative BT methods: maximum change in bounds on two consecutive iterations (default is 1e-3).
 - type of optimization problem solved for BT based on convex envelopes of power flows:
    - 1 - LP with only variable bounds (faster but can produce looser bounds)  (default),
    - 2 - LP with variable bounds and one constraint (slower but can produce tighter bounds).
 - type of optimization problem solved for BT based on convex envelopes of bus current injections:
    - 1 - SOCP with simple cone constraints (fastest) (default),
    - 2 - LP with bounds on variables and one linear constraint (slower),
    - 3 - SOCP+LP (slowest but yields tighter bounds).
 - statistical information to be collected:
    - 0 - nothing,
    - 1 - only general information (default),
    - 2 - general information and detailed information for each method.

### Outputs

#### Test system
Structure in MATPOWER format. If bounds on angle differences were tightened, the updated values are put in the corresponding columns of 'branch' field. 
If bounds on voltage magnitude differences were tightened, extra field 'Vdif' is added to MATPOWER structure. It stores lower and upper bounds on
voltage differences as well as parameters of constraint of type |Vj-slope*Vi|<offset obtained from the properties of the feasible set of the corresponding 
thermal limit constraint. In addition, information related to voltage magnitude differences is stored in columns 22-27 of 'branch' field of the MATPOWER structure. 

#### Statistical information (optional)
Structure containing statistical information about the BT results, e.g. the computation time, percentage of tightened bounds, distance between the bounds etc.

The detailed description of the algorithm's inputs and outputs can be found in function `tighten_bounds.m`.

- - - -

## Related publications
The detailed desciption of the bound tightening algorithm is given in the following publications:
1.  Dmitry Shchetinin, Tomas Tinoco De Rubira, Gabriela Hug
    ["Efficient Bound Tightening Techniques for Convex Relaxations of AC Optimal Power Flow,"](https://ieeexplore.ieee.org/document/8667386) *IEEE
    Transactions on Power Systems*, 2019.  
    DOI: [10.1109/TPWRS.2019.2905232](https://doi.org/10.1109/TPWRS.2019.2905232).

2.  Dmitry Shchetinin, 
    ["Optimization of Power System Operation: Approximations, Relaxations, and Decomposition,"](https://www.research-collection.ethz.ch/handle/20.500.11850/317137) *Doctoral Thesis*, Zurich, ETH Zurich, 2018.
    DOI: [DOI: 10.3929/ethz-b-000317137](https://doi.org/10.3929/ethz-b-000317137).

Please refer to these publications if you use BTOPF in your work.

- - - -

## Example usage
Coming soon

- - - -

## License
BSD 3-clause license.

- - - -

## Authors
Dmitry Shchetinin

- - - -

## Thanks
[Daniel Molzahn](https://molzahn.github.io/) and [Tomas Tinoco De Rubira](https://ttinoco.github.io/) for their valuable suggestions.
