%% load data
dFFT = readOldMriDataOrgData('../mridata/P14/kspace');

d = ifft( ifftshift( dFFT, 3 ), [], 3 );
d = d( :, :, 5:10:end, : );
d = d(:, :, 4, :);

%% convolutions

a = magic(5);

conv_kernel = [2 2 2; 3 0 3; 4 1 2];
%conv_kernel = [2 2; 3 1]
numfits_x = size(a, 2) - size(conv_kernel, 2) + 1;
numfits_y = size(a, 1) - size(conv_kernel, 1) + 1;

out_mx = zeros([numfits_y numfits_x]);
for ix = 1:numfits_x
  for jy = 1:numfits_y
    suba = a(jy:jy+2, ix:ix+2);
    out_mx(jy, ix) = sum(suba.*conv_kernel, 'all');
  end
end

pad_kernel = padData(conv_kernel, size(a));

c = xcorr2(a, conv_kernel);

c1 = conv2(a, flipud(fliplr(conv_kernel)));

%% conv test

fc = flipud(fliplr(conv_kernel));
fcPadded = padData(fc, size(a));

fftfcPadded = fftshift(fft2(ifftshift(fcPadded)));
fftA = fftshift(fft2(ifftshift(a)));
o = ifft2(ifftshift(fftA.*fftfcPadded));
fftshift(o)

test2 = circConv(a, fc)
