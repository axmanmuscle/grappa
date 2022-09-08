function mask = grappa_samplingmask(data_sz, acr_sz, reduction, direction)
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

% we need to make sure that we have the correct center
% if our array is odd size, we choose the center point
% o o o x o o o
% if our array is even size, we choose length / 2 + 1
% o o o x o o 
center_y = ceil((ny+1)/2); 
center_x = ceil((nx+1)/2);

mask = zeros( ny, nx );
switch(direction)
  case 'horiz'
    mask(:,1:reduction:end ) = 1;
  case 'vert'
    mask(1:reduction:end,: ) = 1;
  case 'both'
    mask(1:reduction:end,1:reduction:end ) = 1; % in both directions
  otherwise
    error('Incorrect direction specified. Options: horiz, vert, both');
end
%

mask( center_y-acr_dy:center_y+acr_dy, center_x-acr_dx:center_x+acr_dx ) = 1;
end