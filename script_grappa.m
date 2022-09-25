% now change things around for the new get_kernels
%% load data
dFFT = readOldMriDataOrgData('../mridata/P14/kspace');

d = ifft( ifftshift( dFFT, 3 ), [], 3 );
d = d( :, :, 5:10:end, : );
d = d(:, :, 4, :);
%d = d(:, 5:10:end, :, :);
%d = d(:, 8, :, :);

%% Test Grappa
d1 = squeeze(d);
kernel_sz = [5 5];
acr_sz = [31 31];
reduction = 2;

sm = grappa_samplingmask(size(d1), acr_sz, reduction, 'both');

d1 = bsxfun( @times, d1, sm );

out = grappa(d1, kernel_sz, acr_sz);

r1 = mri_reconSSQ(out);
figure; imshowscale(r1);

r2 = mri_reconSSQ(squeeze(d));
figure; imshowscale(r2);

%% test FFT
d1 = squeeze(d);
kernel_sz = [3 3];
acr_sz = [31 31];
reduction = 2;

sm = grappa_samplingmask(size(d1), acr_sz, reduction, 'horiz');

d1 = bsxfun( @times, d1, sm );

out = grappa_fft(d1, kernel_sz, acr_sz);

%% rotation/ACR testing
d1 = squeeze(d);
%d1 = d1(2:end, 2:end, :);
d2 = rot90(d1, -1);

acr_sz = [15 15];
reduction = 2;

sm = grappa_samplingmask(size(d1), acr_sz, reduction, 'horiz');
d1 = bsxfun( @times, d1, sm );

sm2 = grappa_samplingmask(size(d1), acr_sz, reduction, 'vert');
d2 = bsxfun( @times, d2, sm2 );

kernel_sz = [3 3];
k1 = d1(:, :, 1);
[kers1, karray1] = get_kernels(k1, kernel_sz);
%[os, karrays, kerss] = get_kernel_struct(k1, kernel_sz);

kernel_sz2 = [7 1];
k2 = d2(:, :, 1);
[kers2, karray2] = get_kernels(k2, kernel_sz2);

acr1 = get_acr(d1, acr_sz);
acr2 = get_acr(d2, acr_sz);

kernel = 85;
ka1 = bin_to_array(kernel, kernel_sz);
ka2 = bin_to_array(kernel, kernel_sz2);
[W, A, B] = get_weights(acr1, ka1);
[W2, A2, B2] = get_weights(acr2, ka2);

norm(W - W2)

% fill points i guess
karray_tmp = (karray1 == kernel);
oi = fill_points(d1, karray_tmp, ka1, W);

karray_tmp = (karray2 == kernel);
oi2 = fill_points(d2, karray_tmp, ka2, W2);

o1 = oi(:, :, 1);
o2 = oi2(:, :, 2);
norm(o1 - o2)

out1 = d1 + oi;
out2 = d2 + oi2;

r1 = mri_reconSSQ(out1);
r2 = mri_reconSSQ(out2);

figure; imshowscale(r1);
figure; imshowscale(r2);

%%
%out = grappa(d1, kernel_sz, acr_sz, reduction);
%norm(smash_weights - grappa_weights)

r2 = mri_reconSSQ(squeeze(d));
figure; imshowscale(r2);

r3 = mri_reconSSQ(d1);
figure; imshowscale(r3);

%%  weights testing
d1 = squeeze(d);
%k1 = [1 1 1; 1 0 1; 1 1 1];
%k1 = [1;0;1];
k1 = [1 0 1];

acr_sz = [15 15];
reduction = 2;

acr = get_acr(d1, acr_sz);
weights = get_weights(acr, k1);

%% axial slicing
dFFT = readOldMriDataOrgData('../mridata/P14/kspace');

d = ifft( ifftshift( dFFT, 1 ), [], 1 );
d = d(40, :, :, :);
d1 = squeeze(d);

kernel_sz = [5 5];
acr_sz = [319 31];
reduction = 2;

sm = grappa_samplingmask(size(d1), acr_sz, reduction, 'both');

d1 = bsxfun( @times, d1, sm );

out = grappa(d1, kernel_sz, acr_sz);

r1 = mri_reconSSQ(out);
figure; imshowscale(r1);

r2 = mri_reconSSQ(squeeze(d));
figure; imshowscale(r2);

figure; imshowscale(abs(r1 -r2));

norm(abs(r1 - r2)) / norm(r2)

%% whatever the other slicing is called
dFFT = readOldMriDataOrgData('../mridata/P14/kspace');

d = ifft( ifftshift( dFFT, 2 ), [], 2 );
d = d(:, 40, :, :);
d1 = squeeze(d);

kernel_sz = [5 5];
acr_sz = [31 31];
reduction = 2;

sm = grappa_samplingmask(size(d1), acr_sz, reduction, 'both');

d1 = bsxfun( @times, d1, sm );

out = grappa(d1, kernel_sz, acr_sz);

r1 = mri_reconSSQ(out);
figure; imshowscale(r1);

r2 = mri_reconSSQ(squeeze(d));
figure; imshowscale(r2);

figure; imshowscale(abs(r1 -r2));

norm(abs(r1 - r2)) / norm(r2)