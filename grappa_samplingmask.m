function mask = grappa_samplingmask(data_sz, acr_sz, reduction)
%grappa_samplingmask Creates the desired sampling mask for our grappa alg
%   Detailed explanation goes here
ny = data_sz(1);
nx = data_sz(2);

acr_y = acr_sz(1);
acr_x = acr_sz(2);

if mod(acr_x, 2) ~= 1 || mod(acr_y, 2) ~= 1
  error('ACR Size should be odd in both directions');
end

acr_dy = (acr_y - 1)/2;
acr_dx = (acr_x - 1)/2;

cy = (ny - mod(ny, 2))/2;
cx = (nx - mod(nx, 2))/2;

mask = zeros(data_sz);
mask(:,2:reduction:end-1, :) = 1;
%mask(2:reduction:end-1,:, :) = 0;

mask(cy-acr_dy:cy+acr_dy, cx-acr_dx:cx+acr_dx, :) = 1;
end