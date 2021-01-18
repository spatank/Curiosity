clc; close all; clear;

%% Compute compressibility

data_path = '/Volumes/My Passport/Curiosity/Data/Triangle_Simulations';
files = dir(fullfile(data_path, '0*_*.mat'));

ps = 0.1:0.1:0.9;
iters = 5;

setting = 7;
num_pairs = 100;
compressibility_all = zeros(1, length(ps) * iters);

for i = 1:length(files)
    fprintf('File: %s.\n', files(i).name);
    load(fullfile(data_path, files(i).name));
    fprintf('p = %0.2f, iter = %d.\n', p, iter);
    try
        [S, S_low, clusters, Gs] = rate_distortion_upper_info_new(adj, setting, num_pairs);
        compressibility_all(i) = mean(S(end) - S);
    catch
        compressibility_all(i) = NaN;
    end
end

compressibility = reshape(compressibility_all, length(ps), iters);