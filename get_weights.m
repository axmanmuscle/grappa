function [W, A, B] = get_weights(acr, kernel)
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

% centering
center_y = (kernel_dy - 1)/2;
center_x = (kernel_dx - 1)/2;

% number of sliding kernel fits in the ACR
numfits_x = acr_dx - kernel_dx + 1;
numfits_y = acr_dy - kernel_dy + 1;

% we'll solve B = A * W for W, where B is the collection of points of the
% ACR, A are the corresponding points of the kernel, and W are the weights

% B = [ (numfits_x * numfits_y) x ncoils ]
% W = [ (ncoils * npoints) x ncoils ]
% A = [ (numfits_x * numfits_y) x (ncoils * npoints) ]

% unroll the parallelization to just do one coil at a time
W = zeros([ncoils*npoints ncoils]);
for coilIndx = 1 : ncoils
  A = zeros( [numfits_x * numfits_y ncoils * npoints] );
  B = zeros( [numfits_x * numfits_y 1 ] );
  b_idx = 1;
  for kx = 1:numfits_x
    for ky = 1:numfits_y
      B(b_idx) = acr( ky+center_y, kx+center_x, coilIndx );
      ai = get_points([ky+center_y kx+center_x], acr, kernel);
      A(b_idx, :) = ai;
      b_idx = b_idx + 1;
    end
  end

  w = A \ B;
  W(:, coilIndx) = w;
end

alex = false;
if alex == true
  A = zeros([numfits_x*numfits_y ncoils*npoints]);
  B = zeros([numfits_x*numfits_y ncoils]);
  b_idx = 1;
  for kx = 1:numfits_x
    for ky = 1:numfits_y
      btmp = acr(ky+center_y, kx+center_x, :);
      B(b_idx, :) = reshape(btmp, [ncoils 1]);
      ai = get_points([ky+center_y kx+center_x], acr, kernel);
      A(b_idx, :) = ai;
      b_idx = b_idx + 1;
    end
  end

  W2 = A \ B;
end

