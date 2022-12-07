This repository contains code for handling uncertainty in the projection geometry for computed tomography. It is based on a Bayesian Markov Chain Monte Carlo methodology, and it is connected to the paper 'A Bayesian Approach to CT Reconstruction with Uncertain Geometry'.

For the experiments in the paper, we ran Matlab version 2020a on a linux system.

The below procedures describe how to run the experiments to produce the results of the paper. Slight differences due to different choices of random seed may be observed. The specific resulting data that were used for the paper are provided at
https://doi.org/10.5281/zenodo.7408467
If using these provided data, all steps described below to run experiments can be skipped and just the plotting carried out.

Recreate similar figures (slight deviations due to seed) to the ones in the paper:

Figure 2:

First run the script 'BeadsRecon_Driver.m' to compute and save the reconstructions. Then run the script 'PlotBeads_MAP.m' to plot the reconstructions.

Figure 3-5:

First run the script 'Simple_MCMC_driver.m' for the different configurations (high-dose, low-dose, and short-scan) to compute and save the parameter chains. Then run the script 'ChainPlot.m' to plot the chains.

Figure 6:
Run the script 'Simple_MCMC_driver.m' for all configurations to compute and save the reconstructions. Then run the script 'ReconPlot_AllMethods.m'. 

If you want to play around with other parameters, you can use the script 'Advanced_MCMC_driver.m'.

Dependencies:

astra-toolbox: https://www.astra-toolbox.com/ (v. 1.8.3)

spot-toolbox: https://github.com/mpf/spot

sophiabeads: https://sophilyplum.github.io/sophiabeads-datasets/

SparseBeads data: https://zenodo.org/record/290117#.Y48-4nbMIuU

Nvidia GPU with atleast CUDA 8.0
