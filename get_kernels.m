function [out_struct, karray] = get_kernels(array, kernel_sz)
%get_kernels find all of the kernels we'll need for a given kernel size
% We know in advance how large our kernel is. While most of the points
% needing to be filled in will have the same set of points around them to
% use, we will have edge cases (in most cases, literally at the edge of
% collected k-space). This routine identifies all combinations of points
% which will be used to fill in our reconstruction.
%
% Author: Alex McManus
% *********************
%   Input Parameters:
% *********************
%
%    array:  A 2D (size ny x nx) array of measured k-space. The assumption
%    here is that each coil will have the same undersampling pattern.
%    Uncollected points are assumed to be zero.
%
%    kernel_sz: a 2 element vector containing the dimensions of the kernel.
%    We use MATLAB dimensions, so it's [row column]
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
%    out_struct: a structure that will have numel(out_struct) = # of unique
%    kernels. out_struct(i) will be a matrix of dimensions kernel_sz that
%    will have 1's at collected points and 0's at uncollected points.
%
%    This structure only has one field - ker. This is accessed as
%           out_struct(i).ker
%    which will be the aformentioned kernel matrix.
%    
%    karray: a 2D array of size (ny x nx). This will have 0's at
%    *collected* points. A point with a nonzero value of i will use kernel
%    out_struct(i).ker.

karray = zeros(size(array));

yrows = find(sum(array == 0, 2));
xcols = find(sum(array == 0, 1));

sKerY = kernel_sz(1);
sKerX = kernel_sz(2);
kcell = {};

if mod(sKerY, 2) ~= 1 || mod(sKerX, 2) ~= 1
  error('Kernel dimensions must be odd');
end

ker_dy = (sKerY - 1) / 2;
ker_dx = (sKerX - 1) / 2;

for yi = 1:numel(yrows)
  y = yrows(yi);
  for xi = 1:numel(xcols)
    x = xcols(xi);
    if array(y, x) ~= 0
      continue
    end
    k1 = zeros([sKerY, sKerX]);
    for ky = -ker_dy:ker_dy
      if y + ky < 1 || y + ky > size(array, 1)
        continue
      end

      for kx = -ker_dx:ker_dx
        if x + kx < 1 || x + kx > size(array, 2)
          continue
        end
        if array(y+ky, x+kx) ~= 0
          k1(ky + ker_dy + 1, kx + ker_dx + 1) = 1;
        end

      end
    end
    
    n = numel(kcell);
    found = false;
    for cellidx = 1:n
      if kcell{cellidx} == k1
        found = true;
        karray(y, x) = cellidx;
        break;
      end
    end

    if ~found
      kcell(n+1) = {k1};
      karray(y, x) = n+1;
    end

  end
end

out_struct = struct('ker', kcell);
end