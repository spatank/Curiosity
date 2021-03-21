clc; close all; clear;

%% Preparatory

base_path = '/Volumes/My Passport/Curiosity/';
addpath(genpath(fullfile(base_path, 'Helper')))
addpath(genpath('/Users/sppatankar/Documents/MATLAB/humanStructureFunction'))
data_path = fullfile(base_path, 'v3/Data/reordered_KNOT_preprocessed_Eirene');
files = dir(fullfile(data_path, 'subj_*.mat'));

setting = 7;
num_pairs = 100;

% how large is the largest network?
max_size = 0;
for i = 1:length(files)
    load(fullfile(data_path, files(i).name));
    if size(adj, 1) > max_size
        max_size = size(adj, 1);
    end
end

iters = size(reordered_adjs, 3);

%% Compile metrics of interest for growing knowledge networks
 
% initialize empty matrix of compressibility values
compressibilities = NaN(length(files), iters, max_size);
degrees_of_freedom = NaN(length(files), iters, max_size);

for i = 1:length(files)
    load(fullfile(data_path, files(i).name));
    n = length(adj);
    % adj, reordered_adjs, reordered_nodes, reordered_weighted_Gs available
    for j = 1:iters
        fprintf('Subject %d of %d. Iteration %d of %d.\n', ...
            i, length(files), j, iters)
        G = reordered_weighted_Gs(:, :, j);
        for k = 1:n
            G_filt = zeros(n, n);
            G_filt(1:k, 1:k) = G(1:k, 1:k);
            G_filt(G_filt == 2 * n) = 0; % set 0 weight edges to 0
            G_filt(G_filt > 0) = 1; % binarize
            [components, component_sizes] = conncomp(digraph(G_filt), 'Type', 'Weak');
            idx = component_sizes(components) == max(component_sizes);
            largest_G = full(adjacency(subgraph(digraph(G_filt), idx)));
            [~, num_nodes, num_edges] = density_und(largest_G);
            degrees_of_freedom(i, j, k) = (2 * num_nodes) - num_edges;
            try
                [S, S_low, clusters, Gs] = rate_distortion_upper_info_new(largest_G, setting, num_pairs);
                compressibilities(i, j, k) = mean(S(end) - S);
            catch
                compressibilities(i, j, k) = NaN;
            end
        end
    end
end
% 
% save_string = fullfile(base_path, 'v3/Data/partial_reordered_KNOT_processed.mat');
% save(save_string, 'compressibilities', 'degrees_of_freedom');

%% Preparatory

base_path = '/Volumes/My Passport/Curiosity/';
addpath(genpath(fullfile(base_path, 'Helper')))
addpath(genpath('/Users/sppatankar/Documents/MATLAB/humanStructureFunction'))
data_path = fullfile(base_path, 'v3/Data/reordered_KNOT_preprocessed_Eirene');
files = dir(fullfile(data_path, 'subj_*.mat'));

setting = 7;
num_pairs = 100;

% how large is the largest network?
max_size = 0;
for i = 1:length(files)
    load(fullfile(data_path, files(i).name));
    if size(adj, 1) > max_size
        max_size = size(adj, 1);
    end
end

iters = size(reordered_adjs, 3);

%% Incorporate Betti curves

data_path = fullfile(base_path, 'v3/Data/reordered_KNOT_processed_Eirene');
files = dir(fullfile(data_path, 'subj_*.mat'));

betti_dim_1 = NaN(length(files), iters, max_size);
betti_dim_2 = NaN(length(files), iters, max_size);
betti_dim_3 = NaN(length(files), iters, max_size);

for i = 1:length(files)
    load(fullfile(data_path, files(i).name));
    n = length(adj);
    % adj, barcode_array, betti_curves, n, weighted_adj available
    for j = 1:iters
        fprintf('Subject %d of %d. Iteration %d or %d.\n', ...
            i, length(files), j, iters)
        betti_x = zeros(1, size(betti_curves{j, 1}, 1) - 2);
        betti_dim_1_y = zeros(1, size(betti_curves{j, 1}, 1) - 2);
        betti_dim_2_y = zeros(1, size(betti_curves{j, 2}, 1) - 2);
        betti_dim_3_y = zeros(1, size(betti_curves{j, 3}, 1) - 2);
        for k = 2:size(betti_curves{j, 1}) - 1
            % We ignore 2 entries because of how Eirene processes the networks.
            % The first entry corresponds to the situation when the network has
            % no nodes or edges. The last entry corresponds to the situation
            % when the 0 'edges' are added. 
            betti_x(k - 1) = betti_curves{j, 1}(k, 1);
            betti_dim_1_y(k - 1) = betti_curves{j, 1}(k, 2);
            betti_dim_2_y(k - 1) = betti_curves{j, 2}(k, 2);
            betti_dim_3_y(k - 1) = betti_curves{j, 3}(k, 2);
        end
        betti_dim_1(i, j, betti_x) = betti_dim_1_y;
        betti_dim_2(i, j, betti_x) = betti_dim_2_y;
        betti_dim_3(i, j, betti_x) = betti_dim_3_y;
    end
end

load(fullfile(base_path, 'v3/Data/partial_reordered_KNOT_processed.mat')) ;
save_string = fullfile(base_path, 'v3/Data/all_reordered_KNOT_processed.mat');
save(save_string, 'compressibilities', 'degrees_of_freedom', ...
    'betti_dim_1', 'betti_dim_2', 'betti_dim_3');