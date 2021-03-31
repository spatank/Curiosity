clc; close all; clear;

%% Preparatory

base_path = '/Volumes/My Passport/Curiosity/';
addpath(genpath(fullfile(base_path, 'Helper')))
addpath(genpath('/Users/sppatankar/Documents/MATLAB/humanStructureFunction'))
data_path = fullfile(base_path, 'v5/Data/KNOT/Preprocessed/');

subj_ID = 106;
load(strcat(data_path, 'subj_', string(subj_ID), '_preprocessed.mat'));

setting = 7;
num_pairs = 100;

% num_iters = size(nodes_reordered_weighted, 3); % for null models
num_iters = 5; % for null models

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
    curr_weighted_G(1:n+1:end) = 0;
    curr_G = curr_weighted_G;
    fprintf('Iteration %d of %d.\n', i, num_iters);
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
            fprintf('Rewired: error in stage %d of %d.\n', j, num_iters);
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
    fprintf('Iteration %d of %d.\n', i, num_iters);
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
            fprintf('Latticized: error in stage %d of %d.\n', j, num_iters);
        end
    end
end

%% Plot

% figure;
% hold on
% plot(1:n, DoF, 'Color', [0, 0, 0], ...
%     'LineWidth', 2);
% plot(1:n, mean(DoF_edge_rewired), 'LineWidth', 2, ...
%     'Color', [0.8500, 0.3250, 0.0980]);
% plot(1:n, mean(DoF_latticized), 'LineWidth', 2, ...
%     'Color', [0.9290, 0.6940, 0.1250]);
% hold off
% xlabel('Nodes', 'FontSize', 20);
% ylabel('DoF', 'FontSize', 20);
% legend('Original', 'Edges Rewired', 'Latticized', ...
%     'Location', 'NW');
% title(strcat({'Subj. '}, string(subj_ID), {': DoF'}));
% prettify
% 
% figure;
% hold on
% plot(1:n, C, 'Color', [0, 0, 0], ...
%     'LineWidth', 2);
% plot(1:n, mean(C_edge_rewired), 'LineWidth', 2, ...
%     'Color', [0.8500, 0.3250, 0.0980]);
% plot(1:n, mean(C_latticized), 'LineWidth', 2, ...
%     'Color', [0.9290, 0.6940, 0.1250]);
% hold off
% xlabel('Nodes', 'FontSize', 20);
% ylabel('Compressibility', 'FontSize', 20);
% legend('Original', 'Edges Rewired', 'Latticized', ...
%     'Location', 'NW');
% title(strcat({'Subj. '}, string(subj_ID), {': Compressibility'}));
% prettify
