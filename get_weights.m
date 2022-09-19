function [W, A, B] = get_weights(acr, kernel)
%get_weights solve for weights for given kernel
% This function solves for the prediction weights for a given kernel.
% we'll solve B = A * W for W, where B is the collection of points of the
% ACR, A are the corresponding points of the kernel, and W are the weights
%
% B = [ (numfits_x * numfits_y) x ncoils ]
% W = [ (ncoils * npoints) x ncoils ]
% A = [ (numfits_x * numfits_y) x (ncoils * npoints) ]
%
% Author: Alex McManus
% *********************
%   Input Parameters:
% *********************
%
%    acr:  A 3D (size ny x nx x ncoils) array of fully sampled
%    auto-calibration data.
%
%    kernel: a 2D array describing the current kernel. This is typically
%    output from get_kernels.
%
%    For example, a 3x1 kernel will fit inside a 5x5 ACR as:
%                 x o o o o
%                 x o o o o
%                 x o o o o
%                 o o o o o 
%                 o o o o o 
%
% *********************
%   Output Variables:
% *********************
%
%    W: a 2D array of dimensions [ (ncoils * npoints) ncoils ] containing
%    the weights. The weights for coil i are W(:, i).
%    
%    (OPTIONAL) A: a 2D array of dimensions [ (numfits_x * numfits_y) x ncoils ] 
%    containing the points used for each placement of the kernel within the
%    ACR. Used for testing.
%     
%    (OPTIONAL) B: a 2D array of dimensions [ (numfits_x * numfits_y) (ncoils * npoints) ] 
%    containing each of the center points for the kernel. Used for testing.

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

end

