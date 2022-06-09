function A = Build_Forward_Operator(vol_geom,angles,p,dx,dy,sx,sy,detector_width)
%Function that builds the forward operator A

%INPUT:
%volume_geometry: astra volume geometry object
%angles: projection angles in radians
%p: number of rays
%dx: x-coordinate of detector center
%dy: y-coordinate of detector center
%sx: x-coordinate of source
%sy: y-coordinate of source
%detector_width: Width of detector

%OUTPUT:
%A: forward operator

vectors = zeros(length(angles),6);
vectors(:,1) = sin(angles)*sy + cos(angles)*sx;
vectors(:,2) = -cos(angles)*sy + sin(angles)*sx;
vectors(:,3) = -sin(angles)*dy + dx*cos(angles);
vectors(:,4) = cos(angles)*dy + dx*sin(angles);
vectors(:,5) = cos(angles)*detector_width/p;
vectors(:,6) = sin(angles)*detector_width/p;
    
proj_geom = astra_create_proj_geom('fanflat_vec',p,vectors);
A = opTomo('line_fanflat',proj_geom,vol_geom);
end