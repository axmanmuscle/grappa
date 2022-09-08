function [out_struct, karray, kers] = get_kernel_struct(array, kernel_sz)

% write the inverse function of bin -> arrays

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

%     b = to_dig(k1);
%     %k2 = to_array(b, kernel_sz);
%     %karray(y, x) = b;
% 
%     if isempty(find(kers == b))
%       kers = [kers; b];
%       %b;
%       %k2;
%     end

  end
end
out_struct = struct('ker', kcell);
end

function out = to_dig(array)
  array = reshape(array, 1, []);
  out = bin2dec(num2str(array));
end

function out = to_array(dig, sz)
  if numel(sz) > 2
    error('too many elements in sz');
  end
  bin = dec2bin(dig);
  n = sz(1) * sz(2);
  n1 = n - numel(bin);

  if n1 < 0
    error("doesn't fit");
  end

  tz = repelem('0', n1);
  b1 = strcat(tz, bin);
  out = zeros(sz);

  for i = 1:n
    out(i) = str2num(b1(i));
  end
end