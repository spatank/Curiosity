clc; close all; clear;

%% Specify paths and miscellaneous settings

base_path = '/Volumes/My Passport/Curiosity/';
addpath(genpath(fullfile(base_path, 'Helper')))
addpath(genpath(fullfile(base_path, 'v4/Data')))
data_path = fullfile(base_path, 'v4/Data/KNOT/Raw');
files = dir(fullfile(data_path, 'subj_*.mat'));

%% Build networks

iters = 100; % number of null networks
rewire = 10; % each edge rewired approximately this many times

for i = 1:length(files)
    load(fullfile(data_path, files(i).name))
    % adj, nodes, subj available
    orig_G = double(adj); % some helper functions need double type arguments
    [components, component_sizes] = conncomp(digraph(orig_G), 'Type', 'Weak');
    idx = component_sizes(components) == max(component_sizes);
    G = full(adjacency(subgraph(digraph(orig_G), idx)));
    nodes_kept = nodes(idx, :);
    n = size(G, 1);
    weighted_G = make_weighted_from_order(G, 1:n);
    edges_rewired_weighted = zeros(n, n, iters);
    new_node_order = zeros(n, iters);
    nodes_reordered_weighted = zeros(n, n, iters);
    for j = 1:iters
        edges_rewired_weighted(:, :, j) = ...
            make_weighted_from_order(randmio_und(G, rewire), 1:n);
        new_node_order(:, j) = randperm(n);
        nodes_reordered_weighted(:, :, j) = ...
            make_weighted_from_order(G, new_node_order(:, j));
    end
    save_string = fullfile(base_path, 'v4/Data/KNOT/Preprocessed/', ...
        strcat('subj_', string(subj), '_preprocessed.mat'));
    save(save_string, 'subj', 'nodes', 'orig_G', 'nodes_kept', 'G', 'weighted_G', ...
        'edges_rewired_weighted', 'new_node_order', 'nodes_reordered_weighted');
end