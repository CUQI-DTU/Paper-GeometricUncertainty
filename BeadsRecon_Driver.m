f = fullfile('/zhome','94','f','108663','Desktop','CT-article');
addpath(genpath(f));

%Load the data
beadset = 'B3L1'; % SparseBeads dataset identifier.
setup.beadset = beadset;
toppathname = '/zhome/94/f/108663/Desktop/CT-article/Geometric-Uncertainty-for-X-ray-CT-main/Sophilyplum/SparseBeadsData'; % Parent directory to dataset directory

% Derived parameters
filename = ['SparseBeads_',beadset]; % Name of the dataset e.g. 'SparseBeads_B1L1'
pathname = fullfile(toppathname,filename); % Name of path where the dataset is stored e.g. '/media/somefolder/SparseBeads_B1L1/'.
geom_type = '2D'; % Necessary for loading data. Type can only be '2D' for SparseBeads.

[b,geom] = load_nikon(pathname, filename, geom_type);
b = double(b');
b = reshape(b,size(b,1)*size(b,2),1);

%Set up the geometry
N = 2000;
p = 2000;
theta = 0:0.142857142857143:(2520-1)*0.142857142857143;
detector_width = 400;
angles = geom.angles;
rectSize = 17.4164212609242;

COR = 1;
%cor_vec = 0:0.025:0.3;
%cor_vec = [0,1/10,3/20,0.2,5/20,3/10];
cor_vec = [0,8,11,12,13,16]*(17.4*10^(-3));
dx_vec = 0:1.5:18;

x_recon = zeros(N^2,length(cor_vec));
alpha = 0.12;

%Fista Algorithm parameters
options.maxiters = 500;
options.x0 = zeros(N^2,1);
options.h = 1;
options.epsilon = 10^(-8);
options.iid = 1;
options.nonneg = 1;
options.rho = 10^(-5);
D = speye(N^2);

sy = 121.932688713074;
dy = 1400.204 - sy;

vol_geom = astra_create_vol_geom(N,N,-rectSize,rectSize,-rectSize,rectSize);

if COR == 1
for k=1:length(cor_vec)
    cor = cor_vec(k);
    vectors = zeros(length(angles),6);
    vectors(:,1) = cos(angles)*cor-sin(angles)*sy; 
    vectors(:,2) = sin(angles)*cor+cos(angles)*sy;  
    vectors(:,3) = cos(angles)*cor+sin(angles)*dy;  
    vectors(:,4) = sin(angles)*cor-cos(angles)*dy;  
    vectors(:,5) = cos(angles)*detector_width/p; 
    vectors(:,6) = sin(angles)*detector_width/p;
    
    proj_geom = astra_create_proj_geom('fanflat_vec',p,vectors);
    A = opTomo('cuda',proj_geom,vol_geom);
    normA = normest(A);
    
    x_recon(:,k) = lsqr_gentikh_fixedNorm(A,D,b,zeros(N^2,1),1,alpha,normA,options);
end
else
   for k=1:length(cor_vec)
    dx = dx_vec(k);
    vectors = zeros(length(angles),6);
    vectors(:,1) = cos(angles)*0-sin(angles)*sy; 
    vectors(:,2) = sin(angles)*0+cos(angles)*sy;  
    vectors(:,3) = cos(angles)*dx*0.2+sin(angles)*dy;  
    vectors(:,4) = sin(angles)*dx*0.2-cos(angles)*dy;  
    vectors(:,5) = cos(angles)*detector_width/p; 
    vectors(:,6) = sin(angles)*detector_width/p;
    
    proj_geom = astra_create_proj_geom('fanflat_vec',p,vectors);
    A = opTomo('cuda',proj_geom,vol_geom);
    normA = normest(A);
    
    x_recon(:,k) = lsqr_gentikh_fixedNorm(A,D,b,zeros(N^2,1),1,alpha,normA,options);
   end
end

folder_dir = fullfile('/zhome','94','f','108663','Desktop','CT-article','Results','Real_Data','Beads');

filename = datestr(now, 'dd-mm-yy-HH:MM:SS');

f = fullfile(folder_dir,filename);
save(f,'x_recon','cor_vec','dx_vec')

