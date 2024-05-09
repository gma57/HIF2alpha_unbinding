![HIF2alpha](exit_final.png)
# How does a ligand exit from a buried receptor cavity? Atomistic simulations of unbinding pathways with rigorous kinetics
Marion L. Silvestrini *, Riccardo Solazzo *, Matteo Masetti, Kevin H. Gardner, Lillian T. Chong

## Abstract
Hypoxia-inducible factors (HIFs) are heterodimeric transcription factors that can promote cancer growth. The development of small-molecule drugs that can inhibit the formation of the dimer is therefore a promising route to the treatment of cancer. Here, we focus on the relevant domain of the protein, the HIF2α PAS-B domain, which contains a preformed, buried cavity that binds  artificial small-molecule ligands that allosterically perturb the formation of the HIF heterodimer.  We examine how a representative ligand (THS-017) dissociates and re-enters the buried cavity using atomistic simulations. To enable these simulations, we applied the weighted ensemble path sampling strategy, which can generate continuous pathways with rigorous kinetics (i.e., rate constants) in orders of magnitude less computing time compared to standard simulations. Results reveal a diverse set of pathways for both the ligand unbinding and rebinding processes with estimated rate constants and methyl order parameters that are consistent with experiment. 

### Copyright 

# Weighted Ensemble simulation of the HIF2alpha PAS-B ligand unbinding process

Authors: Riccardo Solazzo, Gessica Adornato

E-mail: riccardo.solazzo@phys.chem.ethz.ch

## Introduction
The purpose of this notebook and the associated github repository is to demonstrate the simulation of the unbinding process of THS-017 from HIF2alpha using the Weighted Ensemble (WE) algorithm in WESTPA 2.0. By executing the provided code, users can replicate the results outlined in the associated publication (DOI: insert DOI).

There are three main steps included in the notebook: Conventional molecular dynamics, weighted ensemble, and analysis. This repository was designed for the user to start at whichever step they want. You can start at the beginning and create your own basis states, or if you want to jump straight to WE, the repo already contains sample basis states. There is also a sample .h5 output file for those who want to see how analysis was conducted without having to run their own WE simulation.

The recommended way to use this notebook is to enter the commands directly into a LINUX terminal. Access to computational resources such as a HPC is also encouraged.

## Table of contents
- [Environment](#environment)
- [Step 1: Generating basis states](#step-1-generating-basis-states)
- [Step 2: Running WE](#step-2-running-we)
    - [Submitting jobs on HPC ](#submitting-jobs-on-hpc)
    - [Simulation monitoring](#simulation-monitoring)
        - [Log file](#log-file)
        - [WEDAP](#wedap)
- [Step 3: Simulation analysis](#step-3-simulation-analysis)
    - [Checking convergence to steady state](#checking-convergence-to-steady-state)

## Environment

First, make sure that all the required packages needed for running the simulations are properly installed.
In this work, we employed the AMBER MD engine, so in order to reproduce our simulation you will need a version of this software ready to run on GPU. Please refer to the AMBER website for instructions on installing AMBER in your particular computing environment.
In addition, the WESTPA package must be installed in your environment using pip:

~~~bash
# Create a conda environment for WESTPA
$ conda create -n westpa --yes python=3.9
# Activate the environment
$ conda activate westpa
# Install WESTPA through pip
$ pip install westpa
~~~

In order to run WESTPA in your computing environment, the environment must be set up by properly configuring the env.sh script. Please edit the included script according to your needs. Assistance with this, including examples for various computing environments, can be found on the WESTPA wiki: https://github.com/westpa/westpa/wiki.

## Step 1: Generating Basis States
Before running a WE simulation we need to generate basis states. To do that, we are going to run conventional MD simulations of the solvated ligand-bound HIF2alpha system. All the files needed for this process are in the modeling_eq_cMD folder. 

First, we use tleap to generate a solvated structure based on the ligand-bound crystal structure using the selected force fields (OPC/ff19sb). To do this, open the terminal in this folder and type the following command:

~~~bash
$ tleap -f tleap.in >tleap.out
~~~
Check the output file to ensure there were no errors in execution. This command should create three files: bound_new.inpcrd, bound_new.prmtop, and bound_new.pdb. Verify that the pdb file has the correct connectivity through a program such as vmd.

If everything looks good, then it is time to run the equilibration. This is best done on an HPC cluster using the run_equil.slurm script:
~~~bash
$ sbatch run_equil.slurm
~~~
This will automatically run through four equilibration steps: minimization, heating, pressurizing, and equilibration. When they are finished, check the output to be sure it ran correctly. There should be 4 new folders in the directory with a .nc, .nfo, .out, and .rst in each.

Next, edit run.slurm to your specifications. prodmanager.sh will use this script to run five 100 ns production runs automatically:
~~~bash
$ bash prodmanager.sh
~~~
The production runs may take a day or two to complete. When they are done, 

## Step 2: Running WE

~~~bash
# Execute the env.sh script to set up the environment variables.
$ ./env.sh
# Execute the init.sh script to initialize the weighted ensemble simulation.
$ ./init.sh
~~~
#### Submitting jobs on HPC

The following command can be used for running the simulation on your local machine:

~~~bash
$ ./run.sh
~~~

In order to run the simulation on a cluster (recommended) a slurm script needs to be created. Once you requested the proper amount of computational resources the simulation can be submitted running the following command:

~~~bash
$ sbatch runwe2.slurm
~~~

### Simulation Monitoring

#### Log file
The easiest way to check the status of the simulation is to read the log file generated during the simulation.

~~~bash
$ tail -n 50 west.log
~~~

#### WEDAP
To have an idea of the amount of configuration space explored so far it can be useful to visualize the probability distribution of your data. To do so we will use wedap a library written by Darian Yang.

In the following cell we will visualize the 2D probability distribution of the single components of the multidimensional progress coordinate, the points are colored according to the negative logarithm of the weights (prior to reweight).

~~~bash
# Install WEDAP in your conda environment
$ pip install wedap
~~~

~~~python
#Import WEDAP to calculate the 2D probability distributions and visualize the results.
import matplotlib.pyplot as plt
import numpy as np
import wedap

wedap.H5_Plot(h5="../../../../HIF/HIF2alpha_unbinding/west.h5", data_type="evolution").plot()
plt.xlabel("Progress Coordinate")
plt.ylabel("WE Iteration")
plt.show()
~~~

It can be useful to plot the first and the second dimensions of the progress coordinate one against the other. In this way you will get an idea of the correlation between two fundamental events for the unbinding mechanism of THS-017: (1) The opening of the buried pocket and (2) the detachment of the ligand from the protein.

~~~python
import matplotlib.pyplot as plt
import numpy as np
import wedap

#wedap.H5_Plot(h5="your_west.h5", data_type="average", Xindex=0, Yindex=1).plot()
wedap.H5_Plot(h5="../../../../HIF/HIF2alpha_unbinding/west.h5", data_type="average", Xindex=0, Yname='IRMSD').plot()
plt.xlabel("Progress Coordinate")
#plt.ylabel("Second Dimension of the progress coordinate")
plt.ylabel("Auxiliary Dataset")
plt.show()
~~~

## Step 3: Simulation Analysis

### Checking convergence to steady state
A WE simulation is converged when it relaxes to a steady state flux. In order to visualize the flux as a function of the molecular time we will need to run the w_ipa command.
You will need to edit the w_ipa settings in the west.cfg file, then you can run the w_ipa command

~~~bash
$ cd ../../HIF/HIF2alpha_unbinding/
$ w_ipa -ra
~~~

The ANALYSIS folder contains the direct.h5 which can be used to visualize the conditional probability fulx.

Finally we can plot the results.

~~~python
import numpy as np
import matplotlib.pyplot as plt
import h5py

direct = h5py.File('../../../../HIF/HIF2alpha_unbinding/ANALYSIS/OVERALL/direct.h5', 'r')

flux_AB=direct['conditional_flux_evolution'][:,0,1]['expected']
flux_A=direct['conditional_flux_evolution'][:,0,1]['expected']/1e-10


plt.plot(flux_A)
plt.xlabel('WE iteration')
plt.ylabel('k$_{off}$ ($s^{-1}$)')
plt.yscale('log')
~~~
