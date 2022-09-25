function k_out = grappa_fft(k_in, kernel_sz, acr_sz)
%grappa an implementation of GRAPPA using the FFT
% Generalized Autocalibrating Partially Parallel Acquisitions (GRAPPA) is a
% parallel imaging reconstruction algorithm for magnetic resonance
% imaging. It reconstructs missing samples as a linear combination of
% nearby points in k-space from *all available coils* by using an
% autocalibration region (ACR), at the center of the image, to solve for
% the interpolation weights given an interpolation kernel size.
%
% This version uses the Fast Fourier Transform to apply the convolution
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

% get auto-calibration region
acr = get_acr(k_in, acr_sz);

% for a given size of kernel, we'll have edge cases depending on the
% sampling pattern
% for each of these edge cases, we'll have to solve for a different set of
% weights
% we assume the undersampling pattern is the same for each coil, so use the
% first coil 
k1 = squeeze(k_in(: ,:, 1));
[kernels, karray] = get_kernels(k1, kernel_sz);
k_out = k_in;

for i = 1:numel(kernels)
  ka = kernels(i).ker;
  Wi = get_weights(acr, ka);

  karray_temp = (karray == i);
  oi = fill_points_fft(k_in, karray_temp, ka, Wi);
  k_out = k_out + oi;
end


end

%%
% okay so the first thing is fill points. for each kernel, iterate through
% each coil. for each coil, take the set of weights and reform them into
% kernel_sz matrices with 0 in the middle. then do the convolution with the
% respective coils, sum them up, then use those values to fill in the
% uncollected points
%
% the next idea will be to use just a single kernel, solve for the one set
% of weights, and then do that for all points and skip all of the get
% weights nonsense (and the get kernels time)
%
% finally, move from padding to just using the specific values of the
% fourier transform using size2fftCoordinates, etc.