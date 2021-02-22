clc; close all; clear;

%% Network Metrics for KNOT

addpath(genpath('/Users/sppatankar/Documents/MATLAB/BCT'))
base_path = '/Volumes/My Passport/Curiosity/';
addpath(genpath(fullfile(base_path, 'Helper')))
data_path = fullfile(base_path, 'v2/Data/KNOT/KNOT_processed_Eirene/');
files = dir(fullfile(data_path, 'subj_*.mat'));

X = zeros(length(files), 2);

for i = 1:length(files)
    load(fullfile(data_path, files(i).name));
    X(i, 1) = charpath(adj);
    X(i, 2) = mean(clustering_coef_bu(adj));
end

k = 3;
idx = kmeans(X, k);

find(X(:, 1) > mean(X(:, 1)))

%% 

save_string = fullfile(base_path, 'v2/Data/KNOT/KNOT_k_means.mat');
save(save_string, 'idx');