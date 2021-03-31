clc; close all; clear;

%% Specify paths and miscellaneous settings

base_path = '/Volumes/My Passport/Curiosity/';
addpath(genpath(fullfile(base_path, 'Helper')))
addpath(genpath(fullfile(base_path, 'v5/Data')))
data_path = fullfile(base_path, 'v5/Data/KNOT/Raw');
files = dir(fullfile(data_path, 'subj_*.mat'));

%% Build networks

iters = 25; % number of null networks
rewire = 50; % each edge rewired approximately this many times

for i = 1:length(files)
    fprintf('Participant %d of %d.\n', i, length(files))
    load(fullfile(data_path, files(i).name))
    % adj, nodes, subj available
    G = double(adj); % some helper functions need double type arguments
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
    save_string = fullfile(base_path, 'v5/Data/KNOT/Preprocessed/', ...
        strcat('subj_', string(subj), '_preprocessed.mat'));
    save(save_string, 'subj', 'nodes', 'G', 'weighted_G', ...
        'edges_rewired_weighted', 'latticized_weighted');
end