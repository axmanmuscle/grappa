function out = fill_points(array, kernel_array, kernel, weights)

[rows, cols] = find(kernel_array ~= 0);
out = zeros(size(array));

n = length(rows);
for i = 1:n
  ri = rows(i);
  ci = cols(i);
  pts = get_points([ri ci], array, kernel);

  fillpt = weights.' * squeeze(pts);
  out(ri, ci, :) = fillpt;
end

end