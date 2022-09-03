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
pts = array(py-kdy:py+kdy, px-kdx:px+kdx,:); 
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