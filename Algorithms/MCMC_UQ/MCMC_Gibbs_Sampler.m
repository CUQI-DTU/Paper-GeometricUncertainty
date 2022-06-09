function res = MCMC_Gibbs_Sampler(setup)
%This function samples from the full posterior (precision variables,
%reconstruction and geometric model parameters) using a Metropolis within
%Gibbs MCMC sampler

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
SOURCE_X = setup.SOURCE_X;                      %Boolean specifying if sx should be estimated                      
SOURCE_Y = setup.SOURCE_Y;                      %Boolean specifying if sy should be estimated
DETECTOR_X = setup.DETECTOR_X;                  %Boolean specifying if dx should be estimated
DETECTOR_Y = setup.DETECTOR_Y;                  %Boolean specifying if dy should be estimated
TILT = setup.TILT;

%TRUE MODEL PARAMETERS
sx_true = setup.sx_true;                        %True value of sx                        
sy_true = setup.sy_true;                        %True value of sy
dx_true = setup.dx_true;                        %True value of dx
dy_true = setup.dy_true;                        %True value of dy
tilt_true = setup.tilt_true;                    %True value of tilt

%INITIAL MODEL PARAMETERS
sx0 = setup.sx0;                                %Initial guess of sx                                  
sy0 = setup.sy0;                                %Initial guess of sy
dx0 = setup.dx0;                                %Initial guess of dx
dy0 = setup.dy0;                                %Initial guess of dy
tilt0 = setup.tilt0;                            %Initial guess of tilt

x0 = setup.x0;                                  %Initial guess of x

%MODEL PARAMETER PRIORS
sx_sigma_prior = setup.sx_sigma_prior;          %std. for sx Gaussian prior
sx_mean_prior = setup.sx_mean_prior;            %mean for sx Gaussian prior
sy_sigma_prior = setup.sy_sigma_prior;          %std. for sy Gaussian prior
sy_mean_prior = setup.sy_mean_prior;            %mean for sy Gaussian prior

dx_sigma_prior = setup.dx_sigma_prior;          %std for dx Gaussian prior
dx_mean_prior = setup.dx_mean_prior;            %mean for dx Gaussian prior
dy_sigma_prior = setup.dy_sigma_prior;          %std for dy Gaussian prior
dy_mean_prior = setup.dy_mean_prior;            %mean for dy Gaussian prior

tilt_mean_prior = setup.tilt_mean_prior;        %mean for tilt gaussian prior
tilt_sigma_prior = setup.tilt_sigma_prior;      %std for tilt gaussian prior

%MODEL PARAMETER Metropolis-Hastings PROPOSALS
sx_sigma_proposal = setup.sx_sigma_proposal;    %std for sx MH proposal
sy_sigma_proposal = setup.sy_sigma_proposal;    %std for sy MH proposal
dx_sigma_proposal = setup.dx_sigma_proposal;    %std for dx MH proposal
dy_sigma_proposal = setup.dy_sigma_proposal;    %std for dy MH proposal
tilt_sigma_proposal = setup.tilt_sigma_proposal;%std for tilt MH proposal

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

if SOURCE_X == 1
    sx_samps = zeros(1,N_samples);      %Array for sx samples
    sx = sx0;                           %Initial sx estimate    
    n_accept_sx = zeros(1,N_samples);                    %Metropolis acceptance rate tracker
else
    sx = sx_true;                       %If we do not estimate the parameter, we set it to the true value
end

if SOURCE_Y == 1
    sy_samps = zeros(1,N_samples);      %Array for sy samples
    sy = sy0;                           %Initial sy estimate
    n_accept_sy = zeros(1,N_samples);                    %Metropolis acceptance rate tracker
else
    sy = sy_true;                       %If we do not estimate the parameter, we set it to the true value
end

%DETECTOR DISTANCE
if DETECTOR_X == 1
    dx_samps = zeros(1,N_samples);      %Array for dx samples
    dx = dx0;                           %Initial dx estimate
    n_accept_dx = zeros(1,N_samples);                    %Metropolis acceptance rate tracker
else
    dx = dx_true;                       %If we do not estimate the parameter, we set it to the true value
end
if DETECTOR_Y == 1
    dy_samps = zeros(1,N_samples);      %Array for dy samples
    dy = dy0;                           %Initial dy estimate
    n_accept_dy = zeros(1,N_samples);                    %Metropolis acceptance rate tracker
else
    dy = dy_true;                       %If we do not estimate the parameter, we set it to the true value
end

if TILT == 1
    tilt_samps = zeros(1,N_samples);
    tilt = tilt0;
    n_accept_tilt = zeros(1,N_samples);
else
    tilt = tilt_true;
end

%Set initial x
x = x0;

%Convert projection angles to radians
angles = theta*pi/180;

%Set volume geometry
vol_geom = astra_create_vol_geom(N,N,-rectSize,rectSize,-rectSize,rectSize);

%Specify initial vectorized geometry
vectors = zeros(length(theta),1);
vectors(:,1) = cos(angles)*sx+sin(angles)*sy;    %First component for source location
vectors(:,2) = sin(angles)*sx-cos(angles)*sy;    %Second component for source location
vectors(:,3) = cos(angles)*dx-sin(angles)*dy;    %First component for detector center
vectors(:,4) = sin(angles)*dx+cos(angles)*dy;    %Second component for detection center
vectors(:,5) = cos(angles + tilt*pi/180)*detector_width/p;     %First component of detector basis
vectors(:,6) = sin(angles + tilt*pi/180)*detector_width/p;     %Second component of detector basis

%Generate forward operator using Astra
proj_geom = astra_create_proj_geom('fanflat_vec',p,vectors);
A = opTomo('line_fanflat',proj_geom,vol_geom);

%Precompute forward projection for efficiency
Ax = A*x;
tic
%Start sampling using Hierarchial Gibbs
for k=1:N_samples
    if mod(k,1)==0
        disp(['Gibbs iteration number: ' num2str(k)])
        disp(['Compute Time for Sample: ' num2str(toc)])
    end
    tic
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
    %METROPOLIS-HASTINGS SAMPLING OF MODEL PARAMETERS
    if SOURCE_X == 1
        count = 0;
        %Sample sx parameter using Metropolis-Hastings 
        for j = 1:N_metro
            %Get candidate sample from proposal
            sx_prop = normrnd(sx,sx_sigma_proposal);
        
            vectors = zeros(length(angles),6);
        
            %Specify vectorized geometry
            vectors(:,1) = cos(angles)*sx_prop+sin(angles)*sy;
            vectors(:,2) = sin(angles)*sx_prop-cos(angles)*sy;
            vectors(:,3) = cos(angles)*dx-sin(angles)*dy;
            vectors(:,4) = sin(angles)*dx+cos(angles)*dy;
            vectors(:,5) = cos(angles + tilt*pi/180)*detector_width/p;
            vectors(:,6) = sin(angles + tilt*pi/180)*detector_width/p;
            
            %Get proposed forward operator
            proj_geom = astra_create_proj_geom('fanflat_vec',p,vectors);
            A_prop = opTomo('line_fanflat',proj_geom,vol_geom);
        
            %Ensure more numerical stability by doing logpdf
            lp_old = -lambda/2*norm(Ax-b,2)^2-1/2*(sx-sx_mean_prior)^2/sx_sigma_prior^2;
            Apropx = A_prop*x;
            lp_prop = -lambda/2*norm(Apropx-b,2)^2-1/2*(sx_prop-sx_mean_prior)^2/sx_sigma_prior^2;
        
            %Check if we accept
            if lp_prop-lp_old>log(rand())
                count = count + 1;
                A = A_prop;
                sx = sx_prop;
                Ax = Apropx;
            end
        end
        if mod(k,100)==0
            disp(['New sx sample: ' num2str(sx)])
            disp(['Number of accepted Metropolis for sx sampling: ' num2str(count)])
        end
        sx_samps(k) = sx;
        n_accept_sx(k) = count;
    end
    if SOURCE_Y == 1
        count = 0;
        %Sample sy parameter using Metropolis-Hastings 
        for j = 1:N_metro
            %Get candidate sample from proposal
            sy_prop = normrnd(sy,sy_sigma_proposal);
        
            vectors = zeros(length(angles),6);
        
            %Specify vectorized geometry
            vectors(:,1) = cos(angles)*sx+sin(angles)*sy_prop;
            vectors(:,2) = sin(angles)*sx-cos(angles)*sy_prop;
            vectors(:,3) = cos(angles)*dx-sin(angles)*dy;
            vectors(:,4) = sin(angles)*dx+cos(angles)*dy;
            vectors(:,5) = cos(angles + tilt*pi/180)*detector_width/p;
            vectors(:,6) = sin(angles + tilt*pi/180)*detector_width/p;
            
            %Get proposed forward operator
            proj_geom = astra_create_proj_geom('fanflat_vec',p,vectors);
            A_prop = opTomo('line_fanflat',proj_geom,vol_geom);
        
            %Ensure more numerical stability by doing logpdf
            lp_old = -lambda/2*norm(Ax-b,2)^2-1/2*(sy-sy_mean_prior)^2/sy_sigma_prior^2;
            Apropx = A_prop*x;
            lp_prop = -lambda/2*norm(Apropx-b,2)^2-1/2*(sy_prop-sy_mean_prior)^2/sy_sigma_prior^2;
        
            %Check if we accept
            if lp_prop-lp_old>log(rand())
                count = count + 1;
                A = A_prop;
                sy = sy_prop;
                Ax = Apropx;
            end
        end
        if mod(k,100)==0
            disp(['New sy sample: ' num2str(sy)])
            disp(['Number of accepted Metropolis for sy sampling: ' num2str(count)])
        end
        sy_samps(k) = sy;
        n_accept_sy(k) = count;
    end
    
    if DETECTOR_X == 1
        count = 0;
        %Sample detector x parameter using Metropolis-Hastings 
        for j = 1:N_metro
            %Get candidate sample from proposal
            dx_prop = normrnd(dx,dx_sigma_proposal);
        
            vectors = zeros(length(angles),6);
        
            %Specify vectorized geometry
            vectors(:,1) = cos(angles)*sx+sin(angles)*sy;
            vectors(:,2) = sin(angles)*sx-cos(angles)*sy;
            vectors(:,3) = cos(angles)*dx_prop-sin(angles)*dy;
            vectors(:,4) = sin(angles)*dx_prop+cos(angles)*dy;
            vectors(:,5) = cos(angles + tilt*pi/180)*detector_width/p;
            vectors(:,6) = sin(angles + tilt*pi/180)*detector_width/p;
            
            %Get proposed forward operator
            proj_geom = astra_create_proj_geom('fanflat_vec',p,vectors);
            A_prop = opTomo('line_fanflat',proj_geom,vol_geom);
        
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
    if DETECTOR_Y == 1
        count = 0;
        %Sample detector y parameter using Metropolis-Hastings 
        for j = 1:N_metro
            %Get candidate sample from proposal
            dy_prop = normrnd(dy,dy_sigma_proposal);
        
            vectors = zeros(length(angles),6);
        
            %Specify vectorized geometry
            vectors(:,1) = cos(angles)*sx+sin(angles)*sy;
            vectors(:,2) = sin(angles)*sx-cos(angles)*sy;
            vectors(:,3) = cos(angles)*dx-sin(angles)*dy_prop;
            vectors(:,4) = sin(angles)*dx+cos(angles)*dy_prop;
            vectors(:,5) = cos(angles + tilt*pi/180)*detector_width/p;
            vectors(:,6) = sin(angles + tilt*pi/180)*detector_width/p;
            
            %Get proposed forward operator
            proj_geom = astra_create_proj_geom('fanflat_vec',p,vectors);
            A_prop = opTomo('line_fanflat',proj_geom,vol_geom);
        
            %Ensure more numerical stability by doing logpdf
            lp_old = -lambda/2*norm(Ax-b,2)^2-1/2*(dy-dy_mean_prior)^2/dy_sigma_prior^2;
            Apropx = A_prop*x;
            lp_prop = -lambda/2*norm(Apropx-b,2)^2-1/2*(dy_prop-dy_mean_prior)^2/dy_sigma_prior^2;
        
            %Check if we accept
            if lp_prop-lp_old>log(rand())
                count = count + 1;
                A = A_prop;
                dy = dy_prop;
                Ax = Apropx;
            end
        end
        if mod(k,100)==0
            disp(['New dy sample: ' num2str(dy)])
            disp(['Number of accepted Metropolis for dy sampling: ' num2str(count)])
        end
        dy_samps(k) = dy;
        n_accept_dy(k) = count;
    end
    if TILT == 1
        count = 0;
        %Sample detector y parameter using Metropolis-Hastings 
        for j = 1:N_metro
            %Get candidate sample from proposal
            tilt_prop = normrnd(tilt,tilt_sigma_proposal);
        
            vectors = zeros(length(angles),6);
        
            %Specify vectorized geometry
            vectors(:,1) = cos(angles)*sx+sin(angles)*sy;
            vectors(:,2) = sin(angles)*sx-cos(angles)*sy;
            vectors(:,3) = cos(angles)*dx-sin(angles)*dy;
            vectors(:,4) = sin(angles)*dx+cos(angles)*dy;
            vectors(:,5) = cos(angles + tilt_prop*pi/180)*detector_width/p;
            vectors(:,6) = sin(angles + tilt_prop*pi/180)*detector_width/p;
            
            %Get proposed forward operator
            proj_geom = astra_create_proj_geom('fanflat_vec',p,vectors);
            A_prop = opTomo('line_fanflat',proj_geom,vol_geom);
        
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
        tilt_samps(k) = tilt;
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

if SOURCE_X == 1
    res.sx_samps = sx_samps;
    res.n_accept_sx = n_accept_sx;
end
if SOURCE_Y == 1
    res.sy_samps = sy_samps;
    res.n_accept_sy = n_accept_sy;
end
if DETECTOR_X == 1
    res.dx_samps = dx_samps;
    res.n_accept_dx = n_accept_dx;
end
if DETECTOR_Y == 1
    res.dy_samps = dy_samps;
    res.n_accept_dy = n_accept_dy;
end
if TILT == 1
    res.tilt_samps = tilt_samps;
    res.n_accept_tilt = n_accept_tilt;
end
res.x_samps = x_samps;
res.setup = setup;
end
