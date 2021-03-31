clc; close all; clear;

%% Specify paths and miscellaneous settings

base_path = '/Volumes/My Passport/Curiosity/';
addpath(genpath('/Users/sppatankar/Documents/MATLAB/BCT'))
addpath(genpath(fullfile(base_path, 'Helper')))
addpath(genpath(fullfile(base_path, 'v5/Data/Wiki/Raw')))

topics = ["molecular_biology"; "geometry"; "optics"; "software_engineering"; ...
    "robotics"; "abstract_algebra"];

%% Build networks

iters = 25; % number of null networks
rewire = 50; % each edge rewired approximately this many times

for i = 1:length(topics)
    fprintf('Topic %d of %d.\n', i, length(topics))
    topic = topics(i);
    load(strcat(topic, '.mat'));
    % adj, nodes, topics, edge_info available
    G = full(adj); % networks are generated in sparse form
    % something in the generation process makes these networks not be binary
    G(G ~= 0) = 1; % force edges to be 1s and non-edges to be 0s
    n = size(G, 1);
    weighted_G = make_weighted_from_order(G, 1:n);
    edges_rewired_weighted = zeros(n, n, iters);
    latticized_weighted = zeros(n, n, iters);
    for j = 1:iters
        edges_rewired_weighted(:, :, j) = ...
            make_weighted_from_order(randmio_und(G, rewire), 1:n);
        [latticized, ~, ~, ~] = latmio_und(G, rewire);
        latticized_weighted(:, :, j) = ...
            make_weighted_from_order(latticized, 1:n);
    end    
    save_string = fullfile(base_path, 'v5/Data/Wiki/Preprocessed/', ...
        strcat(topic, '_preprocessed.mat'));
    save(save_string, 'topic', 'nodes', 'G', 'weighted_G', ...
        'edges_rewired_weighted', 'latticized_weighted');
end