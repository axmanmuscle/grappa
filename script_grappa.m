
dFFT = readOldMriDataOrgData('../mridata/P14/kspace');

d = ifft( ifftshift( dFFT, 3 ), [], 3 );
d = d( :, :, 5:10:end, : );
d = d(:, :, 4, :);
d1 = squeeze(d);
%out = minismash_x(d1);
kernel_sz = [3 3];
acr_sz = [5 5];
reduction = 2;

sm = grappa_samplingmask(size(d1), acr_sz, reduction);

d1 = d1 .* sm;
out = grappa(d1, kernel_sz, acr_sz, reduction);

out_2 = grappa2(d1, kernel_sz, acr_sz);

b = 1;
% to do:
% add in extra reduction sizes
% change size of kernel (?)
% make sure ACR size works for mask and grappa
% set up extra systems to deal with edge points (how do we index)
% documentation

r_truth = mri_reconSSQ(squeeze(d));
r_under = mri_reconSSQ(d1);
r_recon = mri_reconSSQ(out);

figure; title('truth');
imshowscale(r_truth);

figure; title('undersampled');
imshowscale(r_under);

figure; title('reconstructed');
imshowscale(r_recon);

r_diff = abs(r_recon - r_truth);
figure; imshowscale(r_diff);