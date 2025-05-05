function corr_array = fastCorr2(array1, array2, dim, gpuOn)% Synthetic data for testing

% Handle dim input and gpu input

%% Synthetic Data for Testing
% array1 = rand(10000,4096)';
% array2 = rand(10000,4096)';
% array2(:,1) = 3 * array1(:,1) - 2; % r(1,1) should be + 1;
% array2(:,2) = -17 * array1(:,2) + 8; % r(1,2) should be - 1;
% % % rest of r should be random between +1 and -1

%% Compute correlations on n dimension
% Normalize arrays 
array1_norm = (array1-mean(array1))./std(array1);%zscore(array1);
array2_norm = (array2-mean(array2))./std(array2);%zscore(array2);

%% 
corr_array = (1/(size(array1_norm,1)-1)).*(array1_norm'*array2_norm);

end