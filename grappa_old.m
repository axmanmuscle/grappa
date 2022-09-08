function k_out = grappa_old(k_in, kernel_sz, acr_sz, reduction)
%grappa implementation of grappa
%   Detailed explanation goes here

flipped = false;
ny = size(k_in, 1);
nx = size(k_in, 2);
ncoils = size(k_in, 3);

if length(find(sum(k_in(:, :, 1), 1))) < ny
  flipped = true;
  for c =1:ncoils
    k_in(:, :, c) = k_in(:, :, c).';
  end
end
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

S = zeros([ncoils*(kernel_dx - 1)*kernel_dy numfits_x*numfits_y]);
X = zeros([ncoils numfits_x*numfits_y]);
s_idx = 1;
for kx = 1:numfits_x
  for ky = 1:numfits_y
    xtmp = acr(ky+1, kx+1, :);
    X(:, s_idx) = reshape(xtmp, [8 1]);
    stmp = acr([ky ky+kernel_dy-1], kx:kx+kernel_dx-1, :);
    si = reshape(stmp, [], 1);
    S(:, s_idx) = si;
    s_idx = s_idx + 1;
  end
end

% solve for weights
W = X * pinv(S);

% compare / test
ki = logical([1 1 1; 0 0 0; 1 1 1]);
Wi = get_weights(acr, ki);

% apply to each missing point
k_out = k_in;

row_idxs = find(sum(k_out(:, :, 1) == 0, 2) > 0);

for yi = 1:length(row_idxs)
  y = row_idxs(yi);
  containsy = any(acr_yidx == y);
  for x = 2:nx-1
    containsx = any(acr_xidx == x);
    if containsy && containsx
      continue
    end
    stmp = k_out([y-1 y+1], x-1:x+1, :);
    si = reshape(stmp, [], 1);
    xnew = W * si;
    k_out(y, x, :) = xnew;
  end
end

if flipped
  for c = 1:ncoils
    k_out(:, :, c) = k_out(:, :, c).';
  end
end

end