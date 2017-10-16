## biosdp

Implementation of the model revision method for Haemoglobin production model as described in the paper: *Occupation measure methods for modelling and analysis of biological hybrid automata*.

# Requirements

- Matlab 2016a or later
- Mosek v7.1 or later

**NOTE:** If you are using a different version of Mosek, the numerical results may vary with respect to those reported in the paper.

# Installation

- Download as a ZIP file or clone this repository, and open Matlab in the folder `biosdp-master`. 
- Install spotless:

```matlab
>> cd 'spotless-master'
>> spot_install
```

- Add to the working path all the folders in this repository (adjust the full path to your case):

```matlab
>> addpath(genpath('.../biosdp-master'))
```

- Add to the working path your own Mosek installation folder (adjust the full path to your case):

```matlab
>> addpath('.../mosek/7/toolbox/r2013aom/');
```

# Execution

Launch the script `run_iterative_optimal_control`:

```matlab
>> cd 'Haemoglobin Production Implementation'
>> run_iterative_optimal_control
```

Launch the script for plotting the results:

```matlab
>> plot_results('results.mat')
```
