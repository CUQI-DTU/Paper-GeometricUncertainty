This repository contains code for handling uncertainty in the projection geometry for computed tomography. It is based on a Bayesian Markov Chain Monte Carlo methodology,
and it is connected to the paper 'A Bayesian Approach to CT Reconstruction with Uncertain Geometry'.

How to use:
To obtain similar sample chains to the ones presented in the paper you can run the script 'simple_MCMC_driver.m' for the different configurations. You can then use the
scripts 'ChainPlots.m' and 'ReconPlots.m' to create plots of the parameter chains and plots of the reconstructions. If you want to vary other parameters, you can
use the script 'advanced_MCMC_driver.m'.

For the experiments in the paper, we ran Matlab version 2020a on a linux system.

Dependencies:

astra-toolbox: https://www.astra-toolbox.com/ (v. 1.8.3)

spot-toolbox: https://github.com/mpf/spot

sophiabeads: https://sophilyplum.github.io/sophiabeads-datasets/

SparseBeads data: https://zenodo.org/record/290117#.Y48-4nbMIuU

Nvidia GPU with atleast CUDA 8.0
