function W = get_weights(acr, kernel)
% get_weights
% as a design decision, we will take in the acr and a *single*
% kernel

% if we have a 5x5 ACR and a 3x3 kernel, we will *only* shift like:
% x x x o o    o x x x o 
% x x x o o    o x x x o
% x x x o o -> o x x x o
% o o o o o    o o o o o 
% o o o o o    o o o o o
% so at a higher level, just break the ACR down to the specific 3x3xcoils
% array and then apply the kernel to that

% auto-calibration region dimensions
acr_dy = size(acr, 1);
acr_dx = size(acr, 2);
ncoils = size(acr, 3);

% kernel dimensions
kernel_dy = size(kernel, 1);
kernel_dx = size(kernel, 2);
npoints = sum(kernel, 'all');

% number of sliding kernel fits in the ACR
numfits_x = acr_dx - kernel_dx + 1;
numfits_y = acr_dy - kernel_dy + 1;

% we'll solve X = W * S for W, where X is the collection of points of the
% ACR, S are the corresponding points of the kernel, and W are the weights

% X = [ ncoils x (numfits_x * numfits_y) ]
% W = [ ncoils x (ncoils * npoints) ]
% S = [(ncoils * npoints) x (numfits_x * numfits_y)]

S = zeros([ncoils*npoints numfits_x*numfits_y]);
X = zeros([ncoils numfits_x*numfits_y]);
s_idx = 1;
for kx = 1:numfits_x
  for ky = 1:numfits_y
    xtmp = acr(ky+1, kx+1, :);
    X(:, s_idx) = reshape(xtmp, [ncoils 1]);
    si = get_points([ky+1 kx+1], acr, kernel);
    S(:, s_idx) = si;
    s_idx = s_idx + 1;
  end
end

W = X * pinv(S);