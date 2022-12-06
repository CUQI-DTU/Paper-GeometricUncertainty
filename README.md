This repository contains code for handling uncertainty in the projection geometry for computed tomography. It is based on a Bayesian Markov Chain Monte Carlo methodology, and it is connected to the paper 'A Bayesian Approach to CT Reconstruction with Uncertain Geometry'.

For the experiments in the paper, we ran Matlab version 2020a on a linux system.

Recreate similar figures (slight deviations due to seed) to the ones in the paper:

Figure 2:

First run the script 'BeadsRecon_Driver.m' to compute and save the reconstructions. Then run the script 'PlotBeads_MAP.m' to plot the reconstructions.

Figure 3-5:

First run the script 'Simple_MCMC_driver.m' for the different configurations (high-dose, low-dose, and short-scan) to compute and save the parameter chains. Then run the script 'ChainPlot.m' to plot the chains.

Figure 6:
Run the script 'Simple_MCMC_driver.m' for all configurations to compute and save the reconstructions. Then run the script 'ReconPlot_AllMethods.m'.

We also provided results files for the numerical experiments that we ran in the paper. These can be found in the results folder. You can load these and run the plotting scripts outlined in the above to exactly recreate figures 3-6. 

If you want to play around with other parameters, you can use the script 'Advanced_MCMC_driver.m'.

Dependencies:

astra-toolbox: https://www.astra-toolbox.com/ (v. 1.8.3)

spot-toolbox: https://github.com/mpf/spot

sophiabeads: https://sophilyplum.github.io/sophiabeads-datasets/

SparseBeads data: https://zenodo.org/record/290117#.Y48-4nbMIuU

Nvidia GPU with atleast CUDA 8.0
