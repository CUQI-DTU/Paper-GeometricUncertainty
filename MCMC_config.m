function res = MCMC_config(BeadSet,s,k_gibbs,angular_range,path_beads)
%Wrapper function to simulate the configurations in the article
%INPUT:
%BeadSet: strong indicating the beadset to run the simulations
%s: center-of-rotation proposal step size
%k_gibbs: Total number of samples
%angular_range: Integer indicating the final projection angle
%path_beads: Parent directory to SparseBeadsData

%OUTPUT:
%res: struct containing the results.

%Load the data
beadset = BeadSet; % SparseBeads dataset identifier.
setup.beadset = beadset;
toppathname = path_beads;  % Parent directory to dataset directory

% Derived parameters
filename = ['SparseBeads_',beadset]; %Name of the dataset e.g. 'SparseBeads_B1L1'
pathname = fullfile(toppathname,filename); %Name of path where the dataset is stored e.g. '/media/somefolder/SparseBeads_B1L1/'.
geom_type = '2D'; %Necessary for loading data. Type can only be '2D' for SparseBeads.

[b,geom] = load_nikon(pathname, filename, geom_type);
b = b';
nproj = size(b,1); npixel = size(b,2);

%Cut down the data if neccessary
proj_factor = 1; pixel_factor = 1;
proj_cut = nproj/proj_factor;
pixel_cut = npixel/pixel_factor;

if proj_factor~=1 || pixel_factor ~= 1
    b = downsample_beads(b,proj_cut,pixel_cut);
end

%Consider the limited angle problem, where we only have data for the angles
%0 - fin_angle
fin_angle = angular_range;
n_proj = round(2520*(fin_angle/360));
b = b(1:n_proj,:);

setup.b = double(reshape(b,size(b,1)*size(b,2),1));
b = setup.b;
%Set parameters
setup.N = pixel_cut;
setup.p = pixel_cut;
setup.pixelsize  = 0.2*pixel_factor;
setup.detector_width = 400;
setup.rectSize = 17.4164212609242;
n_angles   = proj_cut;
angle_step = 0.142857142857143*proj_factor;
%setup.theta = 0:angle_step:angle_step*(n_angles-1);
%angles = setup.theta*pi/180;
angles = geom.angles(1:proj_factor:n_proj);
setup.theta = angles*180/pi;

%Model parameters
setup.SOURCE_X = 0;                     %Boolean specifying if x-coordinate of SOURCE should be estimated
setup.SOURCE_Y = 0;                     %Boolean specifying if y-coordinate of SOURCE should be estimated
setup.DETECTOR_X = 0;                   %Boolean specifying if x-coordinate of DETECTOR should be estimated    
setup.DETECTOR_Y = 0;                   %Boolean specifying if y-coordinate of DETECTOR should be estimated
setup.TILT = 0;
setup.COR = 1;
setup.virtual_detector = 0;

setup.sy_true = 121.932688713074;
setup.dy_true = 1400.204 - setup.sy_true;
setup.sx_true = 0;
setup.dx_true = 0;
setup.cor_true = 0;
setup.tilt_true = 0;
    
setup.sx0 = 0;                                          %Initial guess for x-coordinate of SOURCE
setup.sy0 = 121.932688713074;                           %Initial guess for y-coordinate of SOURCE
setup.dx0 = 0;                                          %Initial guess for x-coordinate of DETECTOR
setup.dy0 = 1400.204 - 121.932688713074;                %Initial guess for y-coordinate of DETECTOR
setup.cor0 = 0;
setup.tilt0 = 0;

%virtual detector
if setup.virtual_detector == 1
    amp_factor = (setup.sy_true+setup.dy_true)/setup.sy_true;
    setup.detector_width = setup.detector_width/amp_factor;
    setup.dy_true = 0;
    setup.dy0 = 0;
end

%Gibbs sampler parameters
setup.N_samples = k_gibbs;                 %Number of total Gibbs samples
setup.saveRecons = 0;
setup.Nburn = setup.N_samples - 1000;
setup.N_iter = 20;                      %Number of FISTA iterations for x-sample
setup.N_metro = 10;                     %Number of metropolis hastings samples each Gibbs iteration

setup.alpha_delta = 1;                  %Shape parameter for delta prior
setup.beta_delta = 10^(-4);             %Rate parameter for delta prior
setup.alpha_lambda = 1;                 %Shape parameter for lambda prior
setup.beta_lambda = 10^(-4);            %Rate parameter for lambda prior

%Gaussian
setup.sx_sigma_prior = 0.2;             %std for Gaussian Prior for sx
setup.sx_mean_prior = setup.sx0;        %mean for Gaussian Prior for sx

setup.sy_sigma_prior = 5;               %std for Gaussian Prior for sy
setup.sy_mean_prior = setup.sy0;        %mean for Gaussian Prior for sy

setup.dx_sigma_prior = 2;             %std for Gaussian Prior for dx
setup.dx_mean_prior = setup.dx0;        %mean for Gaussian Prior for dx

setup.dy_sigma_prior = 5;               %std for Gaussian Prior for dy
setup.dy_mean_prior = setup.dy0;        %mean for Gasussian Prior for dy

setup.cor_sigma_prior = 0.35;
setup.cor_mean_prior = setup.cor0;

setup.tilt_sigma_prior = 10;
setup.tilt_mean_prior = setup.tilt0;

%Gaussian proposal
setup.sx_sigma_proposal = 0.0005;        %Step size for MH for sx
setup.sy_sigma_proposal = 0.0005;        %Step size for MH for sy
setup.dx_sigma_proposal = 0.0005;        %Step size for MH for dx
setup.dy_sigma_proposal = 0.0005;        %Step size for MH for dy
setup.cor_sigma_proposal = s;
setup.tilt_sigma_proposal = 0.1;

%Initial x-sample parameters
setup.alpha = 0.005;                       %Regularization parameter for initial x-sample
setup.maxiters = 500;                  %Number of iterations for initial x-sample

%Data simulation parameters
setup.noise_level = 0.02;               %Relative measurement Noise level
setup.inverse_factor = 20.74324;        %Factor determining fineness of grid for data simulation

setup.nonneg = 1;                       %Boolean Specifying if nonnegativity
setup.iid = 1;                          %Boolean Specifying if identity matrix for Gaussian Prior

N = setup.N;
p = setup.p;
theta = setup.theta;
rectSize = setup.rectSize;
sx_true = setup.sx_true;
sy_true = setup.sy_true;
dx_true = setup.dx_true;
dy_true = setup.dy_true;
cor_true = setup.cor_true;
tilt_true = setup.tilt_true;
detector_width = setup.detector_width;
sx0 = setup.sx0;
sy0 = setup.sy0;
dx0 = setup.dx0;
dy0 = setup.dy0;
cor0 = setup.cor0;
tilt0 = setup.tilt0;
iid = setup.iid;
nonneg = setup.nonneg;

vectors = zeros(length(theta),6);
vectors(:,1) = cos(angles)*sx0-sin(angles)*sy0; 
vectors(:,2) = sin(angles)*sx0+cos(angles)*sy0;  
vectors(:,3) = cos(angles)*dx0+sin(angles)*dy0;  
vectors(:,4) = sin(angles)*dx0-cos(angles)*dy0;  
vectors(:,5) = cos(angles + tilt0*pi/180)*detector_width/p; 
vectors(:,6) = sin(angles + tilt0*pi/180)*detector_width/p;

%Compute initial projection and volume geometry geometry
vol_geom = astra_create_vol_geom(N,N,-rectSize,rectSize,-rectSize,rectSize);
proj_geom = astra_create_proj_geom('fanflat_vec',p,vectors);
A = opTomo('cuda',proj_geom,vol_geom);
normA = normest(A);
setup.normA = normA;

%Compute initial reconstruction using FISTA
%Gaussian Prior precision matrix
if iid == 1
    D = speye(N^2);
else
    e_vec = ones(N,1);
    L_1D = spdiags([-e_vec 2*e_vec -e_vec], [-1 0 1],N,N);
    L = kron(speye(N),L_1D) + kron(L_1D,speye(N));
    D = chol(L);
end

%Fista Algorithm parameters
options.maxiters = setup.maxiters;
options.x0 = zeros(N^2,1);
options.h = 1;
options.epsilon = 10^(-8);
options.iid = iid;
options.nonneg = nonneg;

alpha = setup.alpha; %Regularization parameter for initial reconstruction

x_MAP = lsqr_gentikh_fixedNorm(A,D,b,zeros(N^2,1),1,alpha,normA,options);
setup.x0 = x_MAP;

%Do Hierarchial Gibbs sampling
res = MCMC_Gibbs_Sampler_Beads(setup);

%Compute mean of precision variables
mean_lambda = mean(res.lambda_samps(res.setup.N_samples-1000:end));
mean_delta = mean(res.delta_samps(res.setup.N_samples-1000:end));

alpha_mean = mean_delta/mean_lambda;
x_MAP_ini = lsqr_gentikh_fixedNorm(A,D,b,zeros(N^2,1),1,alpha_mean,normA,options);
res.X_MAP_ini = x_MAP_ini;

%Compute MAP estimates using pre-processing center of rotation correction
b_reshape = reshape(b,length(angles),p);

dx_com = com_centering(b_reshape,angles,detector_width);
dx_xcorr = xcorr_centering(b_reshape,detector_width);

if setup.COR == 0
    dx_mean = mean(res.dx_samps(res.setup.N_samples - setup.Nburn:end));
    vectors_com = zeros(length(theta),6);
    vectors_com(:,1) = cos(angles)*sx0-sin(angles)*sy0; 
    vectors_com(:,2) = sin(angles)*sx0+cos(angles)*sy0;  
    vectors_com(:,3) = cos(angles)*dx_com+sin(angles)*dy0;  
    vectors_com(:,4) = sin(angles)*dx_com-cos(angles)*dy0;  
    vectors_com(:,5) = cos(angles + tilt0*pi/180)*detector_width/p; 
    vectors_com(:,6) = sin(angles + tilt0*pi/180)*detector_width/p;

    vectors_xcorr = zeros(length(theta),6);
    vectors_xcorr(:,1) = cos(angles)*sx0-sin(angles)*sy0; 
    vectors_xcorr(:,2) = sin(angles)*sx0+cos(angles)*sy0;  
    vectors_xcorr(:,3) = cos(angles)*dx_xcorr+sin(angles)*dy0;  
    vectors_xcorr(:,4) = sin(angles)*dx_xcorr-cos(angles)*dy0;  
    vectors_xcorr(:,5) = cos(angles + tilt0*pi/180)*detector_width/p; 
    vectors_xcorr(:,6) = sin(angles + tilt0*pi/180)*detector_width/p;

    vectors_mcmc = zeros(length(theta),6);
    vectors_mcmc(:,1) = cos(angles)*sx0-sin(angles)*sy0; 
    vectors_mcmc(:,2) = sin(angles)*sx0+cos(angles)*sy0;  
    vectors_mcmc(:,3) = cos(angles)*dx_mean+sin(angles)*dy0;  
    vectors_mcmc(:,4) = sin(angles)*dx_mean-cos(angles)*dy0;  
    vectors_mcmc(:,5) = cos(angles + tilt0*pi/180)*detector_width/p; 
    vectors_mcmc(:,6) = sin(angles + tilt0*pi/180)*detector_width/p;
else
    cor_com = dx_com*sy0/(dy0 + sy0);
    cor_xcorr = dx_xcorr*sy0/(dy0 + sy0);
    cor_mean = mean(res.cor_samps(res.setup.N_samples - setup.Nburn:end));
    vectors_com = zeros(length(theta),6);
    vectors_com(:,1) = cos(angles)*cor_com-sin(angles)*sy0; 
    vectors_com(:,2) = sin(angles)*cor_com+cos(angles)*sy0;  
    vectors_com(:,3) = cos(angles)*cor_com+sin(angles)*dy0;  
    vectors_com(:,4) = sin(angles)*cor_com-cos(angles)*dy0;  
    vectors_com(:,5) = cos(angles + tilt0*pi/180)*detector_width/p; 
    vectors_com(:,6) = sin(angles + tilt0*pi/180)*detector_width/p;

    vectors_xcorr = zeros(length(theta),6);
    vectors_xcorr(:,1) = cos(angles)*cor_xcorr-sin(angles)*sy0; 
    vectors_xcorr(:,2) = sin(angles)*cor_xcorr+cos(angles)*sy0;  
    vectors_xcorr(:,3) = cos(angles)*cor_xcorr+sin(angles)*dy0;  
    vectors_xcorr(:,4) = sin(angles)*cor_xcorr-cos(angles)*dy0;  
    vectors_xcorr(:,5) = cos(angles + tilt0*pi/180)*detector_width/p; 
    vectors_xcorr(:,6) = sin(angles + tilt0*pi/180)*detector_width/p;

    vectors_mcmc = zeros(length(theta),6);
    vectors_mcmc(:,1) = cos(angles)*cor_mean-sin(angles)*sy0; 
    vectors_mcmc(:,2) = sin(angles)*cor_mean+cos(angles)*sy0;  
    vectors_mcmc(:,3) = cos(angles)*cor_mean+sin(angles)*dy0;  
    vectors_mcmc(:,4) = sin(angles)*cor_mean-cos(angles)*dy0;  
    vectors_mcmc(:,5) = cos(angles + tilt0*pi/180)*detector_width/p; 
    vectors_mcmc(:,6) = sin(angles + tilt0*pi/180)*detector_width/p;
end

proj_geom_com = astra_create_proj_geom('fanflat_vec',p,vectors_com);
A_com = opTomo('cuda',proj_geom_com,vol_geom);

proj_geom_xcorr = astra_create_proj_geom('fanflat_vec',p,vectors_xcorr);
A_xcorr = opTomo('cuda',proj_geom_xcorr,vol_geom);

proj_geom_mcmc = astra_create_proj_geom('fanflat_vec',p,vectors_mcmc);
A_mcmc = opTomo('cuda',proj_geom_mcmc,vol_geom);

x_MAP_com = lsqr_gentikh_fixedNorm(A_com,D,b,zeros(N^2,1),1,alpha_mean,normA,options);
x_MAP_xcorr = lsqr_gentikh_fixedNorm(A_xcorr,D,b,zeros(N^2,1),1,alpha_mean,normA,options);
x_MAP_mcmc = lsqr_gentikh_fixedNorm(A_mcmc,D,b,zeros(N^2,1),1,alpha_mean,normA,options);

%Try another alpha-value as well
new_alpha = 0.005;
x_MAP_com_new = lsqr_gentikh_fixedNorm(A_com,D,b,zeros(N^2,1),1,new_alpha,normA,options);
x_MAP_xcorr_new = lsqr_gentikh_fixedNorm(A_xcorr,D,b,zeros(N^2,1),1,new_alpha,normA,options);
x_MAP_mcmc_new = lsqr_gentikh_fixedNorm(A_mcmc,D,b,zeros(N^2,1),1,new_alpha,normA,options);


res.dx_com = dx_com;
res.dx_xcorr = dx_xcorr;
res.X_MAP_com = x_MAP_com;
res.X_MAP_xcorr = x_MAP_xcorr;
res.x_MAP_mcmc = x_MAP_mcmc;
res.x_MAP_com_new = x_MAP_com_new;
res.x_MAP_xcorr_new = x_MAP_xcorr_new;
res.x_MAP_mcmc_new = x_MAP_mcmc_new;

folder_dir = fullfile('/zhome','94','f','108663','Desktop','CT-article','Results','Real_Data','Beads');

filename = datestr(now, 'dd-mm-yy-HH:MM:SS');

f = fullfile(folder_dir,filename);
save(f,'res')