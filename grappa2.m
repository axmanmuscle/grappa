function k_out = grappa2(k_in, kernel_sz, acr_sz)
%grappa2 implementation of grappa
%   Detailed explanation goes here

ny = size(k_in, 1);
nx = size(k_in, 2);
ncoils = size(k_in, 3);

ny = size(k_in, 1);
nx = size(k_in, 2);

% kernel dimensions
kernel_dy = kernel_sz(1);
kernel_dx = kernel_sz(2);

% auto-calibration region dimensions
acr_dy = acr_sz(1);
acr_dx = acr_sz(2);

% number of sliding kernel fits in the ACR
numfits_x = acr_dx - kernel_dx + 1;
numfits_y = acr_dy - kernel_dy + 1;

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

acr_dx = (acr_dx - 1)/2;
acr_dy = (acr_dy - 1)/2;

acr_xidx = [nx/2 - acr_dx : nx/2 + acr_dx];
acr_yidx = [ny/2 - acr_dy : ny/2 + acr_dy];

acr = k_in(ny/2 - acr_dy : ny/2 + acr_dy, nx/2 - acr_dx : nx/2 + acr_dx, :); % assumes the image has even dimensions

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

for i = 1:length(kers)
  ki = kers(i);
  ka = bin_to_array(ki, kernel_sz);
  Wi = get_weights(acr, ka);

  karray_temp = (karray == ki);
  oi = fill_points(k_in, karray_temp, ka, Wi);
  k_out = k_out + oi;
end


end