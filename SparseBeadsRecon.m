clear
clc
close all

%Load Data
disp('Loading Data...')
%%%% Manually modify these variables %%%%
beadset = 'B1L1'; % SparseBeads dataset identifier.
toppathname = 'C:\Users\fhape\Desktop\PhD\Geometric-Uncertainty-for-X-ray-CT-main\Sophilyplum\SparseBeadsData'; % Parent directory to dataset directory e.g. '/media/somefolder/' if this contains e.g. SparseBeads_B1L1/ directory.
iterations = 12; % Number of FISTA iterations to run

% Derived parameters
filename = ['SparseBeads_',beadset]; % Name of the dataset e.g. 'SparseBeads_B1L1'
pathname = fullfile(toppathname,filename); % Name of path where the dataset is stored e.g. '/media/somefolder/SparseBeads_B1L1/'.
geom_type = '2D'; % Necessary for loading data. Type can only be '2D' for SparseBeads.
experiment_name = 'CGLS_Demo'; % For naming purposes... Change to any relevant experiment name or code.

%%%% ------------------------------- %%%%

%setup;    % Setup the mex files.
            % NOTE: You do not need to run this if the mex files already exist.
            %       For more info, type 'help setup' on the command window.

% Load and apply corrections to data
[b,geom] = load_nikon(pathname, filename, geom_type);
%[geom,centre,midpoint] = centre_geom(b, geom);
%%
b = b';
nproj = size(b,1); npixel = size(b,2);
figure; imagesc(reshape(b,nproj,npixel)); colormap gray; colorbar; title('Before Preprocessing')

%Cut down the data
proj_factor = 1; pixel_factor = 4;
proj_cut = nproj/proj_factor;
pixel_cut = npixel/pixel_factor;

b = downsample_beads(b,proj_cut,pixel_cut);

figure; imagesc(reshape(b,proj_cut,pixel_cut)); colormap gray; colorbar; title('After downsample')

%Compute a reconstruction
N = pixel_cut;
p = pixel_cut;
pixelsize  = 0.2*pixel_factor;
detector_width = pixelsize*p;
RectSize = 17.4164212609242;
n_angles   = proj_cut;
angle_step = 0.142857142857143*proj_factor;
theta = 0:angle_step:angle_step*(n_angles-1);
angles = theta*pi/180;

%Compute the centering.
dx_XCORR = cross_corr_centering(b,detector_width);
dx_COM = my_find_center(b,angles,detector_width,0);

sy = 121.932688713074;
dy = 1400.204 - sy;
sx = 0;
dx = 0;

%Compute reconstructions
%Nominel geometry
vectors = zeros(length(theta),1);
vectors(:,1) = cos(angles)*0 + sin(angles)*sy;
vectors(:,2) = sin(angles)*0 - cos(angles)*sy;
vectors(:,3) = cos(angles)*0 - sin(angles)*dy;
vectors(:,4) = sin(angles)*0 + cos(angles)*dy;
vectors(:,5) = cos(angles)*detector_width/p;
vectors(:,6) = sin(angles)*detector_width/p;

%COM
vectors_COM = zeros(length(theta),1);
vectors_COM(:,1) = cos(angles)*sx + sin(angles)*sy;
vectors_COM(:,2) = sin(angles)*sx - cos(angles)*sy;
vectors_COM(:,3) = cos(angles)*dx_XCORR - sin(angles)*dy;
vectors_COM(:,4) = sin(angles)*dx_XCORR + cos(angles)*dy;
vectors_COM(:,5) = cos(angles)*detector_width/p;
vectors_COM(:,6) = sin(angles)*detector_width/p;

%XCORR
vectors_XCORR = zeros(length(theta),1);
vectors_XCORR(:,1) = cos(angles)*sx + sin(angles)*sy;
vectors_XCORR(:,2) = sin(angles)*sx - cos(angles)*sy;
vectors_XCORR(:,3) = cos(angles)*dx_COM - sin(angles)*dy;
vectors_XCORR(:,4) = sin(angles)*dx_COM + cos(angles)*dy;
vectors_XCORR(:,5) = cos(angles)*detector_width/p;
vectors_XCORR(:,6) = sin(angles)*detector_width/p;

%Sophia method
%vectors_sophia = zeros(length(theta),1);
%vectors_sophia(:,1) = cos(angles)*(-midpoint) + sin(angles)*sy;
%vectors_sophia(:,2) = sin(angles)*(-midpoint) - cos(angles)*sy;
%vectors_sophia(:,3) = cos(angles)*(-midpoint) - sin(angles)*dy;
%vectors_sophia(:,4) = sin(angles)*(-midpoint) + cos(angles)*dy;
%vectors_sophia(:,5) = cos(angles)*detector_width/p;
%vectors_sophia(:,6) = sin(angles)*detector_width/p;

%Sophia method number two
%vectors_sophia2 = zeros(length(theta),1);
%vectors_sophia2(:,1) = cos(angles)*centre + sin(angles)*sy;
%vectors_sophia2(:,2) = sin(angles)*centre - cos(angles)*sy;
%vectors_sophia2(:,3) = cos(angles)*centre - sin(angles)*dy;
%vectors_sophia2(:,4) = sin(angles)*centre + cos(angles)*dy;
%vectors_sophia2(:,5) = cos(angles)*detector_width/p;
%vectors_sophia2(:,6) = sin(angles)*detector_width/p;

%Random guess method
%vectors_guess = zeros(length(theta),1);
%vectors_guess(:,1) = cos(angles)*0 + sin(angles)*sy;
%vectors_guess(:,2) = sin(angles)*0 - cos(angles)*sy;
%vectors_guess(:,3) = cos(angles)*0 - sin(angles)*dy;
%vectors_guess(:,4) = sin(angles)*0 + cos(angles)*dy;
%vectors_guess(:,5) = cos(angles)*detector_width/p;
%vectors_guess(:,6) = sin(angles)*detector_width/p;

vol_geom = astra_create_vol_geom(N,N,-RectSize,RectSize,-RectSize,RectSize);
%vol_geom_guess = astra_create_vol_geom(N,N,-RectSize + centre,RectSize + centre,-RectSize + centre,RectSize + centre);

proj_geom_nom = astra_create_proj_geom('fanflat_vec',p,vectors);
A_nom = opTomo('line_fanflat',proj_geom_nom,vol_geom);

proj_geom_COM = astra_create_proj_geom('fanflat_vec',p,vectors_COM);
A_COM = opTomo('line_fanflat',proj_geom_COM,vol_geom);

proj_geom_XCORR = astra_create_proj_geom('fanflat_vec',p,vectors_XCORR);
A_XCORR = opTomo('line_fanflat',proj_geom_XCORR,vol_geom);

%proj_geom_sophia = astra_create_proj_geom('fanflat_vec',p,vectors_sophia);
%A_sophia = opTomo('line_fanflat',proj_geom_sophia,vol_geom);

%proj_geom_sophia2 = astra_create_proj_geom('fanflat_vec',p,vectors_sophia2);
%A_sophia2 = opTomo('line_fanflat',proj_geom_sophia2,vol_geom);

%proj_geom_guess = astra_create_proj_geom('fanflat_vec',p,vectors_guess);
%A_guess = opTomo('line_fanflat',proj_geom_guess,vol_geom_guess);
%e_vec = ones(N,1);
%L_1D = spdiags([-e_vec 2*e_vec -e_vec], [-1 0 1],N,N);
%L = kron(speye(N),L_1D) + kron(L_1D,speye(N));
%D = chol(L);
D = speye(N^2);
%%
options.maxiters = 20;
options.x0 = zeros(N^2,1);
options.h = 1;
options.epsilon = 10^(-8);
options.iid = 1;
options.nonneg = 1;

alpha = 0.005;
b = reshape(b,length(theta)*p,1);

disp('Computing Reconstruction...')
x_MAP_nom = lsqr_gentikh(A_nom,D,double(b),zeros(N^2,1),1,alpha,options);
x_MAP_COM = lsqr_gentikh(A_COM,D,double(b),zeros(N^2,1),1,alpha,options);
x_MAP_XCORR = lsqr_gentikh(A_XCORR,D,double(b),zeros(N^2,1),1,alpha,options);
%x_MAP_sophia = lsqr_gentikh(A_sophia,D,double(b),zeros(N^2,1),1,alpha,options);
%x_MAP_sophia2 = lsqr_gentikh(A_sophia2,D,double(b),zeros(N^2,1),1,alpha,options);
%x_MAP_guess = lsqr_gentikh(A_guess,D,double(b),zeros(N^2,1),1,alpha,options);
%%
figure
imagesc(reshape(x_MAP_nom,N,N))
title('Nominel geometry')
colormap gray
colorbar

figure
imagesc(reshape(x_MAP_COM,N,N))
title('COM correction')
colormap gray
colorbar

figure
imagesc(reshape(x_MAP_XCORR,N,N))
title('XCORR correction')
colormap gray
colorbar

%figure
%imagesc(reshape(x_MAP_sophia,N,N))
%title('sophia correction')
%colormap gray
%colorbar

%figure
%imagesc(reshape(x_MAP_sophia2,N,N))
%title('sophia correction mod')
%colormap gray
%colorbar

%figure
%imagesc(reshape(x_MAP_guess,N,N))
%title('Guess correction')
%colormap gray
%colorbar