function [data_new,theta_new,p_new]=downsample_beads(data,cutproj,cutpixels)
%This function cuts down the data, so that we have nproj projections and
%npixels detector pixels. nproj should be divisible with 2520, and npixels
%should be divisible with 2000

%INPUT:
%   data: The original sinogram
%   cutproj: Number of projections in new sinogram
%   cutpixels: number of detector pixels in new sinogram

%OUTPUT:
%   datanew: The cut down data set
%   theta: new angles in degrees
%   number of detector pixels

[nproj,npixels] = size(data);
AngularStep=0.142857142857143;

factor_proj = nproj/cutproj;
factor_pixels = npixels/cutpixels;
AngularStep_new = AngularStep*factor_proj;

data_new = data(1:factor_proj:end,1:factor_pixels:end);
theta_new = 0:AngularStep_new:360;
p_new = cutpixels;

end
