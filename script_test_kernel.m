%% load data
dFFT = readOldMriDataOrgData('../mridata/P14/kspace');

d = ifft( ifftshift( dFFT, 3 ), [], 3 );
d = d( :, :, 5:10:end, : );
d = d(:, :, 4, :);

%%
ker_szs = {};
ker_szs{1} = [3 3];
ker_szs{2} = [7 7];
ker_szs{3} = [7 1];
ker_szs{4} = [5 5];
ker_szs{5} = [5 1];
ker_szs{6} = [1 7];
ker_szs{7} = [1 5];

check = false;
for i = 1:numel(ker_szs)
  kernel_sz = ker_szs{i};
  d1 = squeeze(d);
  
  acr_sz = [15 15];
  reduction = 2;
  
  sm = grappa_samplingmask(size(d1), acr_sz, reduction, 'horiz');
  d1 = bsxfun( @times, d1, sm );
  
  k1 = d1(:, :, 1);
  [kers1, karray1] = get_kernels(k1, kernel_sz);
  [os, karrays, kerss] = get_kernel_struct(k1, kernel_sz);

  if check
    if numel(kers1) ~= numel(os)
      error('number of kernel mismatch');
    end
    for j = 1:numel(kers1)
      test_array_1 = karray1 == kers1(j);
      tk1 = bin_to_array(kers1(j), kernel_sz);

      test_array_2 = karrays == j;
      tk2 = os(j).ker;

      if tk1 ~= tk2
        error('kernels different');
      end

      if test_array_1 ~= test_array_2
        error('arrays different');
      end
    end
  end
end