function acr = get_acr(k_in, acr_sz)
%get_acr get ACR points from data
% The point of this helper function is to retrieve the fully sampled region
% in the middle of k-space - the auto-calibration region (ACR).
% This takes the undersampled k-space data and the size of the ACR and
% returns the values in the ACR to be used for solving for the weights in
% GRAPPA or other algorithms
%
% Author: Alex McManus
% *********************
%   Input Parameters:
% *********************
%
%    k_in:  A 3D (size ny x nx x ncoils) array of measured k-space. The
%    assumption in GRAPPA is that there is a fully sampled region in the
%    center of the image. The center point convention is as follows:
%
%    if our array is odd size, we choose the center point
%           o o o x o o o
%    if our array is even size, we choose length / 2 + 1
%           o o o x o o 
%
%    acr_sz: A 2 element vector containing the dimensions of the ACR. ACR
%    must be odd sized in both directions.
%
% *********************
%   Output Variables:
% *********************
%
%    acr: A 3D (size ky x kx x ncoils) array of k-space values for the
%    auto-calibration region.

ny = size(k_in, 1);
nx = size(k_in, 2);

acr_dy = acr_sz(1);
acr_dx = acr_sz(2);

if mod(acr_dx, 2) ~= 1 || mod(acr_dy, 2) ~= 1
  error('ACR Size should be odd in both directions');
end

acr_dx = (acr_dx - 1)/2;
acr_dy = (acr_dy - 1)/2;

% we need to make sure that we have the correct center
% if our array is odd size, we choose the center point
% o o o x o o o
% if our array is even size, we choose length / 2 + 1
% o o o x o o 
center_y = ceil((ny+1)/2); 
center_x = ceil((nx+1)/2);

ypts = center_y - acr_dy : center_y + acr_dy;
xpts = center_x - acr_dx : center_x + acr_dx;

acr = k_in(ypts, xpts, :); 

if sum(acr == 0, 'all') > 0
  error('Auto-Calibration Region not fully sampled.');
end

end