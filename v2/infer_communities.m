clc; close all; clear;

base_path = '/Volumes/My Passport/Curiosity/';
addpath(genpath(fullfile(base_path, 'Helper')))
addpath(genpath('/Users/sppatankar/Documents/MATLAB/BCT'));
addpath(genpath('/Users/sppatankar/Documents/MATLAB/WSBM_v1.2'));
data_path = fullfile(base_path, 'v2/Data/KNOT/KNOT_processed_Eirene/');
files = dir(fullfile(data_path, 'subj_*.mat'));

%% Pre-processing

% Load all KNOT time series

load('/Volumes/My Passport/Curiosity/v2/Data/KNOT/all_KNOT_processed.mat')
num_subjects = size(compressibilities, 1);
all_lengths = zeros(1, num_subjects);

% Re-sample compressibility curves

compressibilities = compressibilities(:, 2:end);
betti_dim_1 = betti_dim_1(:, 2:end);
betti_dim_2 = betti_dim_2(:, 2:end);
betti_dim_3 = betti_dim_3(:, 2:end);

for i = 1:num_subjects
    subj_betti_dim_1 = betti_dim_1(i, :);
    idx = find(isnan(subj_betti_dim_1) == 1);
    try
        all_lengths(i) = idx(1) - 1;
    catch
        fprintf('Largest network for subject %d.\n', i);
        all_lengths(i) = length(subj_betti_dim_1);
    end
end

new_length = round(mean(all_lengths));
resampled_comp = zeros(num_subjects, new_length);
resampled_betti_dim_1 = zeros(num_subjects, new_length);
resampled_betti_dim_2 = zeros(num_subjects, new_length);
resampled_betti_dim_3 = zeros(num_subjects, new_length);

for i = 1:num_subjects
    subj_comp = compressibilities(i, 1:all_lengths(i));
    subj_betti_dim_1 = betti_dim_1(i, 1:all_lengths(i));
    subj_betti_dim_2 = betti_dim_2(i, 1:all_lengths(i));
    subj_betti_dim_3 = betti_dim_3(i, 1:all_lengths(i));
    if any(isnan(subj_comp))
        fprintf('NaN value for subject %d.\n', i);
    end
    subj_len = all_lengths(i);
    resampled_comp(i, :) = interp1(1:length(subj_comp), subj_comp, ...
        linspace(1, length(subj_comp), new_length), 'pchip');
    resampled_betti_dim_1(i, :) = interp1(1:length(subj_betti_dim_1), subj_betti_dim_1, ...
        linspace(1, length(subj_comp), new_length), 'pchip');
    resampled_betti_dim_2(i, :) = interp1(1:length(subj_betti_dim_2), subj_betti_dim_2, ...
        linspace(1, length(subj_comp), new_length), 'pchip');
    resampled_betti_dim_3(i, :) = interp1(1:length(subj_betti_dim_3), subj_betti_dim_3, ...
        linspace(1, length(subj_comp), new_length), 'pchip');
end

%% Construct adjacency matrices

adj_comp = zeros(num_subjects, num_subjects);
adj_betti_1 = zeros(num_subjects, num_subjects);

for i = 1:num_subjects
    for j = 1:i-1
        adj_comp(i, j) = corr(resampled_comp(i,:)', resampled_comp(j,:)');
        adj_comp(j, i) = adj_comp(i, j);
        adj_betti_1(i, j) = ...
            corr(resampled_betti_dim_1(i,:)', resampled_betti_dim_1(j,:)');
        adj_betti_1(j, i) = adj_betti_1(i, j);
    end
end

%% Compressibility: collect edge weights

% In the WSBM code missing values are represented by NaNs.
% However, functions of the BCT typically need 0s for missing values.
% Hence, conversions of NaNs to 0s are made here.
adj_comp(adj_comp == 0) = NaN;
% Are drugs trivially co-prescribed with themselves?
adj_comp(1:size(adj_comp,1)+1:end) = NaN; % remove diagonal
edge_list_comp = Adj2Edg(adj_comp);
edge_weights = edge_list_comp(:,3);

figure;
hold on
histogram(edge_weights)
plot([mean(edge_weights); mean(edge_weights)], ...
    repmat(ylim', 1, 1), '-r', 'LineWidth', 2)
hold off
xlabel('Edge Weight', 'FontSize', 15);
ylabel('Frequency', 'FontSize', 15);
title('KNOT: Compressibility', 'FontSize', 15);

%% Compressibility: run WSBM

% input parameters for WSBM
W_Distr = 'Normal';
E_Distr = 'Bernoulli';
num_iters = 5;

log_edge_list = edge_list_comp;
log_edge_list(:,3) = log(edge_list_comp(:,3));

for k = 2:1:10 % sweep over number of communities
    ModelInputs = cell(numel(num_iters),1);
    for iter = 1:num_iters
        ModelInputs{iter} = {k, 'W_Distr', W_Distr, 'E_Distr', E_Distr};
    end
    [Best_Model, Scores, Models] = wsbmLooper(log_edge_list, ModelInputs);
    filename = sprintf('k_%d_communities', k);
    save(fullfile('WSBM_results', 'comp', filename), ...
        'Best_Model', 'Scores', 'Models');
end

%% Compressibility: determine k

path_1 = '/Users/sppatankar/Developer/My Passport/Curiosity/v2/WSBM_results/comp';
k = 2:1:10;
log_evidence = zeros(1,length(k));
errors = zeros(1,length(k));

for idx = 1:length(k)
    curr_k = k(idx);
    path_2 = sprintf('k_%d_communities.mat', curr_k);
    load(fullfile(path_1, path_2)); % loads in Scores and some other variables
    log_evidence(idx) = mean(Scores);
    errors(idx) = std(Scores);
end

num_iters = 5;

figure;
e = errorbar(k, log_evidence, errors/sqrt(num_iters), ...
    '.k', 'MarkerSize', 25, 'LineWidth', 0.1,...
    'MarkerEdgeColor','black','MarkerFaceColor','black');
xlabel('Number of Blocks', 'FontSize', 15);
ylabel('Log Likelihood', 'FontSize', 15);
title('Compressibility: k vs. Log Likelihood', 'FontSize', 15);
prettify;

%% Compressibility: get consensus partition

% Optimal number of communities is k = 4
optimal_k = 4;
path_1 = '/Users/sppatankar/Developer/My Passport/Curiosity/v2/WSBM_results/comp';
path_2 = sprintf('k_%d_communities.mat', optimal_k);
load(fullfile(path_1, path_2)); % loads in Models and some other variables

VI_mat = zeros(length(Models), length(Models));
for i = 1:length(Models)
    model_1 = Models{i, 1}.Para.mu;
    for j = 1:length(Models)
        model_2 = Models{j, 1}.Para.mu;
        VI_mat(i, j) = -varInfo(model_1, model_2);
    end
end
[~, best_model_idx] = max(sum(VI_mat, 2)); % index of most central model
[~, partition] = max(Models{best_model_idx, 1}.Para.mu);

[X, Y, INDSORT] = grid_communities(partition); 
figure;
imagesc(adj_comp(INDSORT, INDSORT)); % plot ordered adjacency matrix
hold on; % hold on to overlay community visualization
plot(X, Y, 'r', 'linewidth', 2);  
hold off;