function out = fill_points_fft(array, kernel_array, kernel, weights)
%fill_points helper function to fill in missing k-space points
% This function fills in the missing k-space data for a given kernel
%
% Author: Alex McManus
% *********************
%   Input Parameters:
% *********************
%
%     array: the full 3D (ny x nx x ncoils) kspace data
%
%     kernel_array: a [ny x nx] array that is nonzero at points that use
%     this interpolation kernel and zeros elsewhere
%
%     kernel: a 2D array of 1's and 0's that corresponds to the current
%     GRAPPA kernel
%
%     weights: the previously solved-for interpolation weights for the
%     given kernel
%
% *********************
%   Output Variables:
% ********************* 
%
%    out: a 3D (ny x nc x ncoils) array of reconstructed k-space data, only
%    nonzero at the points specified in kernel_array

% 
% [rows, cols] = find(kernel_array ~= 0);
% out = zeros(size(array));
% 
% n = length(rows);
% for i = 1:n
%   ri = rows(i);
%   ci = cols(i);
%   pts = get_points([ri ci], array, kernel);
% 
%   fillpt = weights.' * squeeze(pts);
%   out(ri, ci, :) = fillpt;
% end

fdata = fftshift( fft2( ifftshift( array ) ) );

nk = nnz(kernel);
out = zeros(size(array));
for fill_coil = 1:size(array, 3)
  o = zeros(size(kernel_array));
  for coil = 1:size(array, 3)
    ki = kernel;
    ki(find(ki)) = weights(nk*(coil-1)+1:nk*coil, fill_coil);
    
    fk = flipud(fliplr(ki));
    fkPadded = padData(fk, size(array(:, :, coil)));
    
    fftFk = fftshift(fft2(ifftshift(fkPadded)));
    
    oi = ifft2( ifftshift( fftFk .* fdata(:, :, coil)));
    oi = fftshift(oi);
    o = o + oi;
  end
  out(:, :, fill_coil) = o.*kernel_array;
end

end