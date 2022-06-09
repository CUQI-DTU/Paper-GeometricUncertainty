function center = xcorr_centering(sinogram,detector_width)
%This function computes a centering correction for 2D fan-beam data based
%on a simple cross-correlation approach. This is again a simple
%pre-processing method, which is very cheap

%References:
%"Center of Rotation Automatic Measurement for fan-beam CT system based on
%sinogram image features" Min Yang et al. (2013)

%Input:
%sinogram: The recorded sinogram data (dimension is number of projections x
%number of detector pixels)

%output:
%center: The corrected centering parameter

[~,n_rays] = size(sinogram);

%Compute row_sum and flip
row_sum = sum(sinogram,1);
row_sum_flip = fliplr(row_sum);

%Compute cross-correlation
[c,lags] = xcorr(row_sum_flip,row_sum,'normalized');

%Find index of largest cross-correlation
[~,ind] = max(c);

%Take data-fitting approach for finer approximation
c_fit = c(ind-5:ind+5);
lags_fit = lags(ind-5:ind+5);

xx = linspace(lags_fit(1),lags_fit(end),200000);
c_spline = spline(lags_fit,c_fit,xx);

%Find the maximum value and return corresponding lag
[~,center_ind] = max(c_spline);

center_subpixel = xx(center_ind);

%normalize to detectorpixel size
detector_pixel_width = detector_width/n_rays;

center = center_subpixel*detector_pixel_width;
center = center/2;
end



