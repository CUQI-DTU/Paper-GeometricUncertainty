function center = com_centering(sinogram,angles,detector_width)
%This function finds the center offset (x_d offset) for parallel (or full
%angle fan-beam) geometries.

%References
%Azavedo et al, "Calculation of the Rotational Centers in Computed
%Tomography Sinograms" (1990)

%INPUT
%sinogram: measured sinogram
%nAngles: number of projections
%nRays: Number of rays per projection
%angles: projection angles in radians
%detector_width: Width of detector

[nAngles,nRays] = size(sinogram);

h11 = 1/nAngles;
h12 = 0; h13 = 0; h23 = 0;
h22 = 2/nAngles; h33 = 2/nAngles;

p = sinogram;
center = 0; xbar = 0; ybar = 0;

%Compute the coordinates of the detector pixel


for i=1:nAngles
    sum = 0; wsum = 0;
    for j=1:nRays
        sum = sum + p(i,j);
        wsum = wsum++p(i,j)*j;
    end
    d = wsum/sum;
    %angle = delta*i;
    angle = angles(i);
    csum = cos(angle);
    ssum = sin(angle);
    center = center + (d*(h11+h12*csum+h13*ssum));
    xbar = xbar + (d*(h12+h22*csum+h23*ssum));
    ybar = ybar + (d*(h13+h23*csum+h33*ssum));
end

detector_pixel_width = detector_width/nRays;
center = detector_pixel_width*(nRays/2-center)+detector_pixel_width/2;
end

        