clc; close all; clear;

%% Preparatory

base_path = '/Volumes/My Passport/Curiosity/';
addpath(genpath(fullfile(base_path, 'Helper')))
addpath(genpath('/Users/sppatankar/Documents/MATLAB/humanStructureFunction'))
data_path = fullfile(base_path, 'v6/Data/KNOT/Preprocessed/');

subj_ID = 101;
load(strcat(data_path, 'subj_', string(subj_ID), '_preprocessed.mat'));

setting = 7;
num_pairs = 100;

num_iters = size(edges_rewired_weighted, 3); % for null models

n = size(G, 1);
DoF = zeros(1, n);
C = zeros(1, n);

%% Compile metrics of interest for original growing knowledge network

for i = 1:n
    fprintf('Nodes %d of %d.\n', i, n);
    G_filt = G(1:i, 1:i);
    % remove nodes of degree zero
    node_degree = sum(G_filt);
    keep_nodes = find(node_degree ~= 0);
    G_filt = G_filt(keep_nodes, keep_nodes);
    if isempty(G_filt)
        DoF(i) = NaN;
        C(i) = NaN;
        continue
    end
    [~, num_nodes, num_edges] = density_und(G_filt);
    DoF(i) = (2 * num_nodes) - num_edges;
    try
        [S, S_low, clusters, Gs] = ...
            rate_distortion_upper_info(G_filt, setting, num_pairs);
        C(i) = mean(S(end) - S);
    catch
        fprintf('Error caught; size. %d\n', i);
        C(i) = NaN;
    end
end

%% Edge rewired networks

C_edge_rewired = zeros(num_iters, n);
DoF_edge_rewired = zeros(num_iters, n);

for i = 1:num_iters
    curr_weighted_G = edges_rewired_weighted(:, :, i);
    curr_weighted_G(curr_weighted_G == Inf) = 0;
    curr_weighted_G(curr_weighted_G ~= 0) = 1;
    curr_weighted_G(1:n+1:end) = 0; % remove the diagonal
    curr_G = curr_weighted_G;
    symmetry_flag = issymmetric(curr_G);
    fprintf('Iteration %d of %d; symmetry: %d.\n', i, num_iters, symmetry_flag);
    for j = 1:n
        G_filt = curr_G(1:j, 1:j);
        % remove nodes of degree zero
        node_degree = sum(G_filt);
        keep_nodes = find(node_degree ~= 0);
        G_filt = G_filt(keep_nodes, keep_nodes);
        if isempty(G_filt)
            DoF_edge_rewired(i, j) = NaN;
            C_edge_rewired(i, j) = NaN;
            continue
        end
        [~, num_nodes, num_edges] = density_und(G_filt);
        DoF_edge_rewired(i, j) = (2 * num_nodes) - num_edges;
        try
            [S, S_low, clusters, Gs] = ...
                rate_distortion_upper_info(G_filt, setting, num_pairs);
            C_edge_rewired(i, j) = mean(S(end) - S);
        catch
            C_edge_rewired(i, j) = NaN;
            fprintf('Rewired: error in iter %d at stage %d.\n', i, j);
        end
    end
end

%% Latticized networks

C_latticized = zeros(num_iters, n);
DoF_latticized = zeros(num_iters, n);

for i = 1:num_iters
    curr_weighted_G = latticized_weighted(:, :, i);
    curr_weighted_G(curr_weighted_G == Inf) = 0;
    curr_weighted_G(curr_weighted_G ~= 0) = 1;
    curr_weighted_G(1:n+1:end) = 0;
    curr_G = curr_weighted_G;
    symmetry_flag = issymmetric(curr_G);
    fprintf('Iteration %d of %d; symmetry: %d.\n', i, num_iters, symmetry_flag);
    for j = 1:n
        G_filt = curr_G(1:j, 1:j);
        % remove nodes of degree zero
        node_degree = sum(G_filt);
        keep_nodes = find(node_degree ~= 0);
        G_filt = G_filt(keep_nodes, keep_nodes);
        if isempty(G_filt)
            DoF_latticized(i, j) = NaN;
            C_latticized(i, j) = NaN;
            continue
        end
        [~, num_nodes, num_edges] = density_und(G_filt);
        DoF_latticized(i, j) = (2 * num_nodes) - num_edges;
        try
            [S, S_low, clusters, Gs] = ...
                rate_distortion_upper_info(G_filt, setting, num_pairs);
            C_latticized(i, j) = mean(S(end) - S);
        catch
            C_latticized(i, j) = NaN;
            fprintf('Latticized: error in iter %d at stage %d.\n', i, j);
        end
    end
end

%% Variables of interest

clearvars -except DoF C DoF_edge_rewired DoF_latticized C_edge_rewired ...
    C_latticized n subj_ID
