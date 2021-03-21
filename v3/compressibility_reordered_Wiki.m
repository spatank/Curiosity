clc; close all; clear;

%% Compressibility for Growing Knowledge Networks

base_path = '/Volumes/My Passport/Curiosity/';
addpath(genpath(fullfile(base_path, 'Helper')))
addpath(genpath(fullfile(base_path, 'Data')))
orig_data_path = fullfile(base_path, 'v2/Data/Wiki/Wiki_processed_Eirene/');
null_data_path = fullfile(base_path, 'v3/Data/reordered_Wiki_preprocessed_Eirene/');

topic = 'geometry';

null_data = load(fullfile(null_data_path, strcat('reordered_', topic, '.mat')));

setting = 7;
num_pairs = 100;

%% Collect metrics of interest for original network

orig_data = load(fullfile(orig_data_path, strcat(topic, '.mat')));
n = length(orig_data.adj);
compressibilities = NaN(1, n);
degrees_of_freedom = NaN(1, n);
G = orig_data.weighted_adj;

for k = 1:n
    fprintf('%d nodes out of %d.\n', k, n);
    G_filt = zeros(n, n);
    G_filt(1:k, 1:k) = G(1:k, 1:k);
    G_filt(G_filt == 2 * n) = 0; % set 0 weight edges to 0
    G_filt(G_filt > 0) = 1; % binarize
    [components, component_sizes] = conncomp(digraph(G_filt), 'Type', 'Weak');
    idx = component_sizes(components) == max(component_sizes);
    largest_G = full(adjacency(subgraph(digraph(G_filt), idx)));
    [~, num_nodes, num_edges] = density_und(largest_G);
    degrees_of_freedom(k) = (2 * num_nodes) - num_edges;
    try
        [S, S_low, clusters, Gs] = rate_distortion_upper_info_new(largest_G, setting, num_pairs);
        compressibilities(k) = mean(S(end) - S);
    catch
        compressibilities(k) = NaN;
    end
end

betti_dim_1 = NaN(1, n);
betti_dim_2 = NaN(1, n);
betti_dim_3 = NaN(1, n);

% adj, barcode_array, betti_curves, n, weighted_adj available
betti_x = zeros(1, size(orig_data.betti_curves{1, 1}, 1) - 2);
betti_dim_1_y = zeros(1, size(orig_data.betti_curves{1, 1}, 1) - 2);
betti_dim_2_y = zeros(1, size(orig_data.betti_curves{2, 1}, 1) - 2);
betti_dim_3_y = zeros(1, size(orig_data.betti_curves{3, 1}, 1) - 2);
for k = 2:size(orig_data.betti_curves{1, 1}) - 1
    % We ignore 2 entries because of how Eirene processes the networks.
    % The first entry corresponds to the situation when the network has
    % no nodes or edges. The last entry corresponds to the situation
    % when the 0 'edges' are added. 
    betti_x(k - 1) = orig_data.betti_curves{1, 1}(k, 1);
    betti_dim_1_y(k - 1) = orig_data.betti_curves{1, 1}(k, 2);
    betti_dim_2_y(k - 1) = orig_data.betti_curves{2, 1}(k, 2);
    betti_dim_3_y(k - 1) = orig_data.betti_curves{3, 1}(k, 2);
end
betti_dim_1(betti_x) = betti_dim_1_y;
betti_dim_2(betti_x) = betti_dim_2_y;
betti_dim_3(betti_x) = betti_dim_3_y;

save_string = fullfile(base_path, 'v3/Data/', ...
    strcat('all_', topic, '_processed.mat'));
save(save_string, 'compressibilities', 'degrees_of_freedom', ...
    'betti_dim_1', 'betti_dim_2', 'betti_dim_3');

%% Collect metrics of interest for reordered networks

n = length(null_data.adj);
iters = size(null_data.reordered_adjs, 3);

% initialize empty matrix of compressibility values
compressibilities = NaN(iters, n);
degrees_of_freedom = NaN(iters, n);

for j = 1:iters
    G = null_data.reordered_weighted_Gs(:, :, j);
    for k = 1:n
        fprintf('%d nodes out of %d, in iteration %d.\n', k, n, j);
        G_filt = zeros(n, n);
        G_filt(1:k, 1:k) = G(1:k, 1:k);
        G_filt(G_filt == 2 * n) = 0; % set 0 weight edges to 0
        G_filt(G_filt > 0) = 1; % binarize
        [components, component_sizes] = conncomp(digraph(G_filt), 'Type', 'Weak');
        idx = component_sizes(components) == max(component_sizes);
        largest_G = full(adjacency(subgraph(digraph(G_filt), idx)));
        [~, num_nodes, num_edges] = density_und(largest_G);
        degrees_of_freedom(j, k) = (2 * num_nodes) - num_edges;
        try
            [S, S_low, clusters, Gs] = rate_distortion_upper_info_new(largest_G, setting, num_pairs);
            compressibilities(j, k) = mean(S(end) - S);
        catch
            compressibilities(j, k) = NaN;
        end
    end
end

data_path = fullfile(base_path, 'v3/Data/reordered_Wiki_preprocessed_Eirene');
null_data_Betti = load(fullfile(data_path, ...
    strcat('partial_processed_reordered_', topic, '.mat')));

betti_dim_1 = NaN(iters, n);
betti_dim_2 = NaN(iters, n);
betti_dim_3 = NaN(iters, n);

% adj, barcode_array, betti_curves, n, weighted_adj available
for j = 1:iters
    fprintf('Iteration %d or %d.\n', j, iters)
    betti_x = zeros(1, size(null_data_Betti.betti_curves{j, 1}, 1) - 2);
    betti_dim_1_y = zeros(1, size(null_data_Betti.betti_curves{j, 1}, 1) - 2);
    betti_dim_2_y = zeros(1, size(null_data_Betti.betti_curves{j, 2}, 1) - 2);
    betti_dim_3_y = zeros(1, size(null_data_Betti.betti_curves{j, 3}, 1) - 2);
    for k = 2:size(null_data_Betti.betti_curves{j, 1}) - 1
        % We ignore 2 entries because of how Eirene processes the networks.
        % The first entry corresponds to the situation when the network has
        % no nodes or edges. The last entry corresponds to the situation
        % when the 0 'edges' are added. 
        betti_x(k - 1) = null_data_Betti.betti_curves{j, 1}(k, 1);
        betti_dim_1_y(k - 1) = null_data_Betti.betti_curves{j, 1}(k, 2);
        betti_dim_2_y(k - 1) = null_data_Betti.betti_curves{j, 2}(k, 2);
        betti_dim_3_y(k - 1) = null_data_Betti.betti_curves{j, 3}(k, 2);
    end
    betti_dim_1(j, betti_x) = betti_dim_1_y;
    betti_dim_2(j, betti_x) = betti_dim_2_y;
    betti_dim_3(j, betti_x) = betti_dim_3_y;
end

save_string = fullfile(base_path, 'v3/Data/', ...
    strcat('all_reordered_', topic, '_processed.mat'));
save(save_string, 'compressibilities', 'degrees_of_freedom', ...
    'betti_dim_1', 'betti_dim_2', 'betti_dim_3');