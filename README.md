This repository contains code for handling uncertainty in the projection geometry for computed tomography. It is based on a Bayesian Markov Chain Monte Carlo methodology, and it is connected to the paper 'A Bayesian Approach to CT Reconstruction with Uncertain Geometry'.

For the experiments in the paper, we ran Matlab version 2020a on a linux system.

How to use:

To obtain figure 2 in the paper, first run the script 'BeadsRecon_Driver.m' to compute the reconstructions, and then run the script 'PlotBeads_MAP.m'.

To obtain similar sample chains to the ones presented in the paper you can run the script 'Simple_MCMC_driver.m' for the different configurations. You can use the
script 'ChainPlots.m' to create chain plots similar to figures 3,4,5 in the paper. When results for all configurations have been obtained, you can use 'ReconPlot_AllMethods.m' to create figures of the reconstructions similar to figure 6 in the paper.

If you want to play around with other parameters, you can use the script 'Advanced_MCMC_driver.m'.

Dependencies:

astra-toolbox: https://www.astra-toolbox.com/ (v. 1.8.3)

spot-toolbox: https://github.com/mpf/spot

sophiabeads: https://sophilyplum.github.io/sophiabeads-datasets/

SparseBeads data: https://zenodo.org/record/290117#.Y48-4nbMIuU

Nvidia GPU with atleast CUDA 8.0
