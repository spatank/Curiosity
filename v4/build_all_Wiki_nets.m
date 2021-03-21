clc; close all; clear;

%% Specify paths and miscellaneous settings

base_path = '/Volumes/My Passport/Curiosity/';
addpath(genpath(fullfile(base_path, 'Helper')))
addpath(genpath(fullfile(base_path, 'v4/Data/Wiki/Raw')))
data_path = fullfile(base_path, 'v2/Data/Wiki/Wiki_processed_Eirene/');

topics = ["molecular_biology"; "geometry"; "optics"; "software_engineering"; ...
    "robotics"; "abstract_algebra"];

%% Build networks

iters = 100; % number of null networks
rewire = 10; % each edge rewired approximately this many times

for i = 1:length(topics)
    fprintf('Topic %d of %d.\n', i, length(topics))
    topic = topics(i);
    load(strcat(topic, '.mat'));
    % adj, nodes, topics, edge_info available
    orig_G = full(adj); % networks are generated in sparse form
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
    save_string = fullfile(base_path, 'v4/Data/Wiki/Preprocessed/', ...
        strcat(topic, '_preprocessed.mat'));
    save(save_string, 'topic', 'nodes', 'orig_G', 'nodes_kept', 'G', 'weighted_G', ...
        'edges_rewired_weighted', 'new_node_order', 'nodes_reordered_weighted');
end