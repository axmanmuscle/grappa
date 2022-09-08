function [k_out, weights] = grappa(k_in, kernel_sz, acr_sz)
%grappa an implementation of GRAPPA
% Generalized Autocalibrating Partially Parallel Acquisitions (GRAPPA) is a
% parallel imaging reconstruction algorithm for magnetic resonance
% imaging. It reconstructs missing samples as a linear combination of
% nearby points in k-space from *all available coils* by using an
% autocalibration region (ACR), at the center of the image, to solve for
% the interpolation weights given an interpolation kernel size.
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
%    kernel_sz: a 2 element vector containing the dimensions of the kernel.
%    We use MATLAB dimensions, so it's [row column]
%    For example, a 3x1 kernel will fit inside a 5x5 ACR as:
%                 x o o o o
%                 x o o o o
%                 x o o o o
%                 o o o o o 
%                 o o o o o 
%    acr_sz: A 2 element vector containing the dimensions of the ACR. ACR
%    must be odd sized in both directions.
%
% *********************
%   Output Variables:
% *********************
%
%    k_out: A 3D (size ny x nx x ncoils) array of k-space values with the
%    reconstructed values filled in.
%     
%    weights: **DEBUG** a set of interpolation weights
%


% kernel dimensions
kernel_dy = kernel_sz(1);
kernel_dx = kernel_sz(2);

% auto-calibration region dimensions
acr_dy = acr_sz(1);
acr_dx = acr_sz(2);

% size checking
if mod(kernel_dx, 2) ~= 1 || mod(kernel_dy, 2) ~= 1
  error('Kernel size should be odd in both directions');
end

if mod(acr_dx, 2) ~= 1 || mod(acr_dy, 2) ~= 1
  error('ACR Size should be odd in both directions');
end

if kernel_dx > acr_dx
  error('Kernel is larger than ACR in kx-dimension');
end

if kernel_dy > acr_dy
  error('Kernel is larger than ACR in ky-dimension');
end

%acr_dx = (acr_dx - 1)/2;
%acr_dy = (acr_dy - 1)/2;

%acr = k_in(ny/2 - acr_dy : ny/2 + acr_dy, nx/2 - acr_dx : nx/2 + acr_dx, :); % assumes the image has even dimensions
acr = get_acr(k_in, acr_sz);

if sum(acr == 0, 'all') > 0
  error('Auto-Calibration Region not fully sampled.');
end

% between the reduction factor, the ACR size and the kernel size, we need
% to figure out how to fill out the matrix we'll use to solve for our
% weights

% there will be issues here - how do we handle edge cases, e.g.
% now use the kernel to slide over the ACR

% we'll solve X = W * S for W, where X is the collection of points of the
% ACR, S are the corresponding points of the kernel, and W are the weights

% X = [ ncoils x (numfits_x * numfits_y) ]
% W = [ ncoils x (ncoils * kernel_dx * kernel_dy) ]
% S = [(ncoils * kernel_dx * kernel_dy) x (numfits_x * numfits_y)]

k1 = squeeze(k_in(: ,:, 1));
[kers, karray] = get_kernels(k1, kernel_sz);
k_out = k_in;

weights = 0;
disp(length(kers))
for i = 1:length(kers)
%for i = 2
  ki = kers(i);
  ka = bin_to_array(ki, kernel_sz);
  Wi = get_weights(acr, ka);

  if i == 1
    weights = Wi;
  end

  karray_temp = (karray == ki);
  oi = fill_points(k_in, karray_temp, ka, Wi);
%   for coil = 1:ncoils
%     wii = Wi(:, coil);
%     oi = fill_points_onecoil(k_in, karray_temp, ka, wii, coil);
%     k_out(:, :, coil) = k_out(:, :, coil) + oi;
%   end
  k_out = k_out + oi;
end


end