function out = get_points(pt_idx, array, kernel)
%get_points helper function to retrieve correct points
% this function is used both to set up the weights and get the appropriate
% points to interpolate with when filling in k-space
%
% Author: Alex McManus
% *********************
%   Input Parameters:
% *********************
%
%     pt_idx: A 2 element vector containing the [row column] of a specific
%     point. Coil index does not matter here since for solving for the
%     interpolation weights and filling in the missing values, we need the
%     values from *every* coil, so we pull from every coil all the time
%     regardless of which coil we're *currently* solving for.
%
%     array: the full 3D (ny x nx x ncoils) kspace data
%
%     kernel: a 2D array of 1's and 0's that corresponds to the current
%     GRAPPA kernel
%
% *********************
%   Output Variables:
% ********************* 
%
%    out: the set of points for the input pt_idx over all the coils

k1 = logical(kernel);

kdy = (size(kernel, 1) - 1) / 2;
kdx = (size(kernel, 2) - 1) / 2;

py = pt_idx(1);
px = pt_idx(2);

% need to add in error checking for the kernel here
dist_ym = py - kdy;
dist_yp = py + kdy;
dist_xm = px - kdx;
dist_xp = px + kdx;

if dist_ym < 1
  % at the top of the array
  ypts = 1:py+kdy;
  cf = 1 - dist_ym; % correction factor
  k1 = k1(cf+1:end, :);
elseif dist_yp > size(array, 1)
  % at the bottom
  ypts = py-kdy:size(array, 1);
  cf = dist_yp - size(array, 1);
  k1 = k1(1:end-cf, :);
else
  ypts = py-kdy:py+kdy;
end
if dist_xm < 1
  % on the left edge
  xpts = 1:px+kdx;
  cf = 1 - dist_xm;
  k1 = k1(:, cf+1:end);
elseif dist_xp > size(array, 2)
  % on the right edge
  xpts = px-kdx:size(array, 2);
  cf = dist_xp - size(array, 2);
  k1 = k1(:, 1:end-cf);
else
  xpts = px-kdx:px+kdx;
end

pts = array(ypts, xpts, :);
k2 = repmat(k1, [1, 1, size(array, 3)]);

out = pts(k2);