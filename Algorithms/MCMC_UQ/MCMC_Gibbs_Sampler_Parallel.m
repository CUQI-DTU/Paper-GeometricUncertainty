function res = MCMC_Gibbs_Sampler_Parallel(setup)
%This function samples from the full posterior (precision variables,
%reconstruction and geometric model parameters) using a Metropolis within
%Gibbs MCMC sampler for parallel-beam geometry

%INPUT:
%setup: struct including the following fields

%Projection parameters
N = setup.N;                                    %Reconstruction Grid Size (N x N)
theta = setup.theta;                            %Projection angles (degrees)
p = setup.p;                                    %Number of detector pixels
detector_width = setup.detector_width;          %Width of detector
rectSize = setup.rectSize;                      %Size of object domain (rectSize x rectSize)

b = setup.b;                                    %sinogram data

N_samples = setup.N_samples;                    %Number of Gibbs samples
N_iter = setup.N_iter;                          %Number of FISTA iterations for x-sample
N_metro = setup.N_metro;                        %Number of Metropolis-Hastings iterations each Gibbs iteration

%Gamma conjugate parameters
alpha_delta = setup.alpha_delta;                %Shape parameter for delta (x-prior parameter) prior
beta_delta = setup.beta_delta;                  %Rate parameter for delta (x-prior parameter) prior
alpha_lambda = setup.alpha_lambda;              %Shape parameter for lambda (noise precision) prior
beta_lambda = setup.beta_lambda;                %Rate parameter for lambda (noise precision) prior

%ESTIMATION VARIABLES
DETECTOR_X = setup.DETECTOR_X;                  %Boolean specifying if dx should be estimated
TILT = setup.TILT;

%TRUE MODEL PARAMETERS
dx_true = setup.dx_true;                        %True value of dx
tilt_true = setup.tilt_true;

%INITIAL MODEL PARAMETERS                                 
dx0 = setup.dx0;                                %Initial guess of dx
tilt0 = setup.tilt0;

x0 = setup.x0;                                  %Initial guess of x

%MODEL PARAMETER PRIORS
dx_sigma_prior = setup.dx_sigma_prior;          %std for dx Gaussian prior
dx_mean_prior = setup.dx_mean_prior;            %mean for dx Gaussian prior

tilt_sigma_prior = setup.tilt_sigma_prior;
tilt_mean_prior = setup.tilt_mean_prior;

%MODEL PARAMETER Metropolis-Hastings PROPOSALS
dx_sigma_proposal = setup.dx_sigma_proposal;    %std for dx MH proposal
tilt_sigma_proposal = setup.tilt_sigma_proposal;

%FISTA ALGORITHM DETAILS
nonneg = setup.nonneg;                          %Boolean specifying if nonnegativity constraints
iid =    setup.iid;                             %Specify if identity matrix for Gaussian precision matrix (alt. is derivative)

M = length(theta)*p;

%Set up FISTA algorithm options
%Precision Matrix
if iid == 1
    %Regularization precision matrix is the identity
    D = speye(N^2);
    L = D;
else  
    %Regularization precision matrix is finite difference approximation of
    %laplace operator with homogenous boundary conditions
    e_vec = ones(N,1);
    L_1D = spdiags([-e_vec 2*e_vec -e_vec], [-1 0 1],N,N);
    L = kron(speye(N),L_1D) + kron(L_1D,speye(N));
    D = chol(L);
end
options.nonneg = nonneg;     %Set non-negativity options
options.maxiters = N_iter;   %Set maximum number of iterations
options.h = 1;               %Set discretization parameter
options.epsilon = 10^(-8);   %Heuristic tolerance parameter
options.iid = iid;           %Specify if identity precision matrix

%Preallocate arrays for samples
x_samps =       zeros(N^2,N_samples);
delta_samps =   zeros(1,N_samples);
lambda_samps =  zeros(1,N_samples);

%Preallocate model parameters arrays and set initial model parameters
if DETECTOR_X == 1
    dx_samps = zeros(1,N_samples);      %Array for dx samples
    dx = dx0;                           %Initial dx estimate
    n_accept_dx = zeros(1,N_samples);   %Array for dx acceptance rate
else
    dx = dx_true;                       %If we do not estimate the parameter, we set it to the true value
end

if TILT == 1
    tilt_samps = zeros(1,N_samples);      %Array for dx samples
    tilt = tilt0;                           %Initial dx estimate
    n_accept_tilt = zeros(1,N_samples);   %Array for dx acceptance rate
else
    tilt = tilt_true;                       %If we do not estimate the parameter, we set it to the true value
end

%Set initial x
x = x0;

%Convert projection angles to radians
angles = theta*pi/180;

%Set volume geometry
vol_geom = astra_create_vol_geom(N,N,-rectSize,rectSize,-rectSize,rectSize);

%Specify initial vectorized geometry
vectors = zeros(length(theta),1);
vectors(:,1) = -sin(angles);                     %First component for source location
vectors(:,2) = cos(angles);                      %Second component for source location
vectors(:,3) = cos(angles)*dx;                   %First component for detector center
vectors(:,4) = sin(angles)*dx;                   %Second component for detection center
vectors(:,5) = cos(angles + tilt*pi/180)*detector_width/p;     %First component of detector basis
vectors(:,6) = sin(angles + tilt*pi/180)*detector_width/p;     %Second component of detector basis

%Generate forward operator using Astra
proj_geom = astra_create_proj_geom('parallel_vec',p,vectors);
A = opTomo('line',proj_geom,vol_geom);

%Precompute forward projection for efficiency
Ax = A*x;

%Start sampling using Hierarchial Gibbs
for k=1:N_samples
    if mod(k,100) == 0
        disp(['Gibbs iteration number: ' num2str(k)])
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Sample lambda and delta from conjugate Gamma distributions
    if nonneg == 1
        N_nonzero = nnz(x);
        lambda = gamrnd(M/2+alpha_lambda,1/(1/2*norm(Ax-b)^2+beta_lambda));
        delta = gamrnd(N_nonzero/2+alpha_delta,1/(1/2*x'*L*x+beta_delta));
    else
        lambda = gamrnd(M/2+alpha_lambda,1/(1/2*norm(Ax-b)^2+beta_lambda));
        delta = gamrnd(N^2/2+alpha_delta,1/(1/2*x'*L*x+beta_delta));
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %JOINT METROPOLIS-HASTINGS SAMPLING OF MODEL PARAMETERS
    count = 0;
    if DETECTOR_X == 1
        %Sample detector x parameter using Metropolis-Hastings 
        for j = 1:N_metro
            %Get candidate sample from proposal
            dx_prop = normrnd(dx,dx_sigma_proposal);
        
            vectors = zeros(length(angles),6);
        
            %Specify vectorized geometry
            vectors(:,1) = -sin(angles);
            vectors(:,2) = cos(angles);
            vectors(:,3) = cos(angles)*dx_prop;
            vectors(:,4) = sin(angles)*dx_prop;
            vectors(:,5) = cos(angles + tilt*pi/180)*detector_width/p;
            vectors(:,6) = sin(angles + tilt*pi/180)*detector_width/p;
            
            %Get proposed forward operator
            proj_geom = astra_create_proj_geom('parallel_vec',p,vectors);
            A_prop = opTomo('line',proj_geom,vol_geom);
        
            %Ensure more numerical stability by doing logpdf
            lp_old = -lambda/2*norm(Ax-b,2)^2-1/2*(dx-dx_mean_prior)^2/dx_sigma_prior^2;
            Apropx = A_prop*x;
            lp_prop = -lambda/2*norm(Apropx-b,2)^2-1/2*(dx_prop-dx_mean_prior)^2/dx_sigma_prior^2;
        
            %Check if we accept
            if lp_prop-lp_old>log(rand())
                count = count + 1;
                A = A_prop;
                dx = dx_prop;
                Ax = Apropx;
            end
        end
        
        if mod(k,100)==0
            disp(['New dx sample: ' num2str(dx)])
            disp(['Number of accepted Metropolis for dx sampling: ' num2str(count)])
        end
        dx_samps(k) = dx;
        n_accept_dx(k) = count;
    end
    
    if TILT == 1
        %Sample detector x parameter using Metropolis-Hastings 
        for j = 1:N_metro
            %Get candidate sample from proposal
            tilt_prop = normrnd(tilt,tilt_sigma_proposal);
        
            vectors = zeros(length(angles),6);
        
            %Specify vectorized geometry
            vectors(:,1) = -sin(angles);
            vectors(:,2) = cos(angles);
            vectors(:,3) = cos(angles)*dx;
            vectors(:,4) = sin(angles)*dx;
            vectors(:,5) = cos(angles + tilt_prop*pi/180)*detector_width/p;
            vectors(:,6) = sin(angles + tilt_prop*pi/180)*detector_width/p;
            
            %Get proposed forward operator
            proj_geom = astra_create_proj_geom('parallel_vec',p,vectors);
            A_prop = opTomo('line',proj_geom,vol_geom);
        
            %Ensure more numerical stability by doing logpdf
            lp_old = -lambda/2*norm(Ax-b,2)^2-1/2*(tilt-tilt_mean_prior)^2/tilt_sigma_prior^2;
            Apropx = A_prop*x;
            lp_prop = -lambda/2*norm(Apropx-b,2)^2-1/2*(tilt_prop-tilt_mean_prior)^2/tilt_sigma_prior^2;
        
            %Check if we accept
            if lp_prop-lp_old>log(rand())
                count = count + 1;
                A = A_prop;
                tilt = tilt_prop;
                Ax = Apropx;
            end
        end
        
        if mod(k,100)==0
            disp(['New tilt sample: ' num2str(tilt)])
            disp(['Number of accepted Metropolis for tilt sampling: ' num2str(count)])
        end
        tilt_samps(k) = dx;
        n_accept_tilt(k) = count;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Sample approximate x from multivariate normal distribution using FISTA
    
    %Get i.i.d normal samples
    e_N = randn(N^2,1);
    e_M = randn(M,1);
    
    %Perturb sinogram
    b_tilde = b + lambda^(-1/2)*e_M;
    u = delta^(-1/2)*e_N;
    
    %Run FISTA starting at previous sample for N_iter iterations
    options.x0 = x;
    x = lsqr_gentikh(A,D,b_tilde,u,lambda,delta,options);
    Ax = A*x;
    %save samples    
    delta_samps(k) = delta;
    lambda_samps(k) = lambda;
    x_samps(:,k) = x;
end

%Save results
res.delta_samps = delta_samps;
res.lambda_samps = lambda_samps;

if DETECTOR_X == 1
    res.dx_samps = dx_samps;
    res.n_accept_dx = n_accept_dx;
end
if TILT == 1
    res.tilt_samps = tilt_samps;
    res.n_accept_tilt = n_accept_tilt;
end

res.x_samps = x_samps;
res.setup = setup;
end