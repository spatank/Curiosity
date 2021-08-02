clc; close all; clear;

%% Specify paths and miscellaneous settings

base_path = '/Volumes/My Passport/Curiosity/';
addpath(genpath('/Users/sppatankar/Documents/MATLAB/BCT'))
addpath(genpath(fullfile(base_path, 'Helper')))
addpath(genpath(fullfile(base_path, 'v8/Data/Wiki/Raw')))

topics = ["molecular_biology"; "geometry"; "optics"; "software_engineering"; ...
    "abstract_algebra"];

%% Build networks

iters = 25; % number of null networks
rewire = 50; % each edge rewired approximately this many times

for i = 1:length(topics)
    fprintf('Topic %d of %d.\n', i, length(topics))
    topic = topics(i);
    load(strcat(topic, '.mat'));
    % adj, nodes, topics, edge_info available
    G = double(adj); % some helper functions need double type arguments
    G(logical(eye(size(G)))) = 0; % set diagonal entries to 0
    n = size(G, 1);
    weighted_G = make_weighted_from_order(G, 1:n);
    edges_rewired_weighted = zeros(n, n, iters);
    latticized_weighted = zeros(n, n, iters);
    for j = 1:iters
        G_rewired = randmio_und(G, rewire);
        % sometimes the rewirings are asymmetric
        while ~issymmetric(G_rewired)
            G_rewired = randmio_und(G, rewire);
        end
        edges_rewired_weighted(:, :, j) = ...
            make_weighted_from_order(G_rewired, 1:n);
        [G_latticized, ~, ~, ~] = latmio_und(G, rewire);
        % sometimes the rewirings are asymmetric
        while ~issymmetric(G_latticized)
            [G_latticized, ~, ~, ~] = latmio_und(G, rewire);
        end
        latticized_weighted(:, :, j) = ...
            make_weighted_from_order(G_latticized, 1:n);
    end    
    save_string = fullfile(base_path, 'v8/Data/Wiki/Preprocessed/', ...
        strcat(topic, '_preprocessed.mat'));
    save(save_string, 'topic', 'nodes', 'G', 'weighted_G', ...
        'edges_rewired_weighted', 'latticized_weighted');
end