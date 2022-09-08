%% Testing grappa components
a = zeros(15);
a(1:2:end, 1:2:end) = 1;
kernel_sz = [5 5];
[kers, karray] = get_kernels(a, kernel_sz);
b = repmat(a, [1 1 3]);
for i = 1:length(kers)
  ki = kers(i);
  kai = karray == ki;

  o = fill_points(a, kai, bin_to_array(ki, kernel_sz), ones(3));
end

%% Testing weights
testA = [17     4    24     6     1    13     3    12    13    10    29    33    24    27    34    32    21    21;
    24     6     1    13     8    20    13    10    29    33     2     9    34    32    21    21     4     2;
     1    13     8    20    15    22    29    33     2     9    31    21    21    21     4     2     2    23;
    23    10     5    12     7    19    15    22     5    37    35     9     8    15    19    33    20     4;
     5    12     7    19    14    21     5    37    35     9     8    40    19    33    20     4    35    25;
     7    19    14    21    16     3    35     9     8    40    19    38    20     4    35    25    32    49;
     4    11     6    18    13    25    12    13    10    11    33     7    27    27    32    21    21     3;
     6    18    13    25    20     2    10    11    33     7     9    15    32    21    21     3     2     7;
    13    25    20     2    22     9    33     7     9    15    21     5    21     3     2     7    23    17];
testX = [5 5 19; 7 35 20; 14 8 35; 6 10 32; 13 33 21; 20 9 2; 12 37 33; 19 9 4; 21 40 25];

testACR = zeros([5 5 3]);
testACR(:, :, 1) = [17    24     1     8    15;
                    23     5     7    14    16;
                     4     6    13    20    22;
                    10    12    19    21     3;
                    11    18    25     2     9];
testACR(:, :, 2) = [3    13    29     2    31;
                   15     5    35     8    19;
                   12    10    33     9    21;
                   22    37     9    40    38;
                   13    11     7    15     5];
testACR(:, :, 3) = [24    34    21     4     2;
                     8    19    20    35    32;
                    27    32    21     2    23;
                    15    33     4    25    49;
                    27    21     3     7    17];
testKernel = [1 1 1; 0 0 0 ; 1 1 1];

testWeights = pinv(testA) * testX;

[fnWeights, fnA, fnB] = get_weights(testACR, testKernel);
if norm(testWeights - fnWeights) > 1e-12
  error("Error in weights");
end
o = get_points([2 2], testACR, testKernel);
tp = squeeze(testACR(2, 2, :));
tp2 = testWeights.' * o;
if norm(tp2 - tp) > 1e-14
  error("Error in filling points");
end