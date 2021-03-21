clc; close all; clear;

base_path = '/Volumes/My Passport/Curiosity/';
addpath(genpath(fullfile(base_path, 'Helper')))
addpath(genpath(fullfile(base_path, 'v3/Data/')))
load('prop_prob_v3.mat') % model
% loads in all_weighted_Gs, bettiCurve, compressibilities


%% Compile metrics of interest for growing knowledge networks

iters = size(all_weighted_Gs, 3);
n = size(all_weighted_Gs, 1);

setting = 7;
num_pairs = 100;

degrees_of_freedom = NaN(iters, n);

for i = 1:iters
    fprintf('Iteration %d of %d\n', i, iters);
    G = all_weighted_Gs(:, :, i);
    for k = 1:n
        G_filt = zeros(n, n);
        G_filt(1:k, 1:k) = G(1:k, 1:k);
        G_filt(G_filt == 2 * n) = 0; % set 0 weight edges to 0
        G_filt(G_filt > 0) = 1; % binarize
        [components, component_sizes] = conncomp(digraph(G_filt), 'Type', 'Weak');
        idx = component_sizes(components) == max(component_sizes);
        largest_G = full(adjacency(subgraph(digraph(G_filt), idx)));
        [~, num_nodes, num_edges] = density_und(largest_G);
        degrees_of_freedom(i, k) = (2 * num_nodes) - num_edges;        
    end
end

betti_dim_1 = bettiCurve(:, 2:end, 1);
betti_dim_2 = bettiCurve(:, 2:end, 2);
betti_dim_3 = bettiCurve(:, 2:end, 3);

save_string = fullfile(base_path, 'v3/Data/', 'osc_prob_processed.mat');
save(save_string, 'compressibilities', 'degrees_of_freedom', ...
    'betti_dim_1', 'betti_dim_2', 'betti_dim_3');
