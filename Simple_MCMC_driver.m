%Driver script to obtain sample chains similar to the ones shown in the
%article. There may be some slight differences due to the seed

%All the units in the software are in physical units of millimeter instead of pixels, which is why
%the numerical values for the center-of-rotation proposal step size is different

%Set your parent directory for the SparseBeads data
path_beads = '/zhome/94/f/108663/Desktop/CT-article/Geometric-Uncertainty-for-X-ray-CT-main/Sophilyplum/SparseBeadsData';

%Configurations from the paper

%high-dose configuration
BeadSet = 'B3L1';
s = 2*10^(-5);
angular_range = 360;
k_gibbs = 6000;

%low-dose configuration
%BeadSet = 'Low_B3L1';
%s = 5*10^(-5);
%angular_range = 360;
%k_gibbs = 6000;

%short-scan configuration
%BeadSet = 'B3L1';
%s = 5*10^(-5);
%angular_range = 210;
%k_gibbs = 10000;

%Perform the experiments (initialization, MCMC sampling, COM/XCORR
%correction)
res = MCMC_config(BeadSet,s,k_gibbs,angular_range,path_beads);

%Save the results
filename = 'B3L1_results'; %Filename for the results file
save(filename,'res')