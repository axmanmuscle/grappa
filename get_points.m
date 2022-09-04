function out = get_points(pt_idx, array, kernel)
%get_points helper function to retrieve correct points
% this function is used both to set up the weights and get the appropriate
% points to interpolate with when filling in k-space

k1 = logical(kernel);

kdy = (size(kernel, 1) - 1) / 2;
kdx = (size(kernel, 2) - 1) / 2;

py = pt_idx(1);
px = pt_idx(2);

% need a try/catch here for the boundary points
% points on the boundary still get a [kernel size] kernel, so if y == 1
% then this will throw an exception

% figure out if we're boundary or not
left = false;
right = false;
top = false;
bottom = false;

% need to add in error checking for the kernel here
if py == 1
  top = true;
  ypts = py:py+kdy;
  k1 = k1(1+kdy:end, :);
elseif py == size(array, 1)
  bottom = true;
  ypts = py-kdy:py;
  k1 = k1(1:1+kdy, :);
else
  ypts = py-kdy:py+kdy;
end
if px == 1
  left = true;
  xpts = px:px+kdy;
  k1 = k1(:, 1+kdy:end);
elseif px == size(array, 2)
  right = true;
  xpts = px-kdx:px;
  k1 = k1(:, 1:1+kdy);
else
  xpts = px-kdx:px+kdx;
end

pts = array(ypts, xpts, :);
k2 = repmat(k1, [1, 1, size(array, 3)]);

out = pts(k2);

%%% may not have to do this
% if we have a 5x5 ACR and a 3x3 kernel, we will *only* shift like:
% x x x o o    o x x x o 
% x x x o o    o x x x o
% x x x o o -> o x x x o
% o o o o o    o o o o o 
% o o o o o    o o o o o
% so at a higher level, just break the ACR down to the specific 3x3xcoils
% array and then apply the kernel to that