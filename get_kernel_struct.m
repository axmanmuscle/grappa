function [out_struct, karray, kers] = get_kernel_struct(array, kernel_sz)

kers = [];
karray = zeros(size(array));

yrows = find(sum(array == 0, 2));
xcols = find(sum(array == 0, 1));

ker_y = kernel_sz(1);
ker_x = kernel_sz(2);
kcell = {};

if mod(ker_y, 2) ~= 1 || mod(ker_x, 2) ~= 1
  error('Kernel dimensions must be odd');
end

ker_dy = (ker_y - 1) / 2;
ker_dx = (ker_x - 1) / 2;

for yi = 1:numel(yrows)
  y = yrows(yi);
  for xi = 1:numel(xcols)
    x = xcols(xi);
    if array(y, x) ~= 0
      continue
    end
    k1 = zeros([ker_y, ker_x]);

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
      %ctmp = ;
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