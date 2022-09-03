function out = bin_to_array(dig, sz)
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