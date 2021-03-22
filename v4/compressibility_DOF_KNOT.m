clc; close all; clear;

%% Preparatory

base_path = '/Volumes/My Passport/Curiosity/';
addpath(genpath(fullfile(base_path, 'Helper')))
addpath(genpath('/Users/sppatankar/Documents/MATLAB/humanStructureFunction'))
data_path = fullfile(base_path, 'v4/Data/KNOT/Preprocessed/');

subj_ID = 106;
load(strcat(data_path, 'subj_', string(subj_ID), '_preprocessed.mat'));

setting = 7;
num_pairs = 100;

%% Compile metrics of interest for original growing knowledge network

n = size(G, 1);
DoF = zeros(1, n);
C = zeros(1, n);

for i = 1:n
    fprintf('Nodes %d of %d.\n', i, n);
    G_filt = G(1:i, 1:i);
    [~, num_nodes, num_edges] = density_und(G_filt);
    DoF(i) = (2 * num_nodes) - num_edges;
    [components, component_sizes] = conncomp(digraph(G_filt), 'Type', 'Weak');
    num_components = length(unique(components));
    compressibilities_filt = zeros(1, num_components);
    for k = 1:num_components
        curr_component = components == k;
        component_subgraph = full(adjacency(subgraph(digraph(G_filt), ...
            curr_component)));
        try
            [S, S_low, clusters, Gs] = ...
                rate_distortion_upper_info_new(component_subgraph, setting, num_pairs);
            compressibilities_filt(k) = mean(S(end) - S);
        catch
            compressibilities_filt(k) = NaN;
        end
    end
    C(i) = sum(compressibilities_filt, 'omitnan');
end

%% Compile metrics of interest for null models of growing knowledge network

% Node re-ordered networks

% num_iters = size(nodes_reordered_weighted, 3);
num_iters = 10;
C_node_reordered = zeros(num_iters, n);
DoF_node_reordered = zeros(num_iters, n);

for i = 1:num_iters
    curr_weighted_G = nodes_reordered_weighted(:, :, i);
    curr_weighted_G(curr_weighted_G == Inf) = 0;
    curr_weighted_G(curr_weighted_G ~= 0) = 1;
    curr_weighted_G(1:n+1:end) = 0;
    curr_G = curr_weighted_G;
    fprintf('Iteration %d of %d.\n', i, 10);
    curr_node_order = new_node_order(:, i);
    for j = 1:n
        G_filt = curr_G(1:j, 1:j);
        [~, num_nodes, num_edges] = density_und(G_filt);
        DoF_node_reordered(i, j) = (2 * num_nodes) - num_edges;
        [components, component_sizes] = conncomp(digraph(G_filt), 'Type', 'Weak');
        num_components = length(unique(components));
        compressibilities_filt = zeros(1, num_components);
        for k = 1:num_components
            curr_component = components == k;
            component_subgraph = full(adjacency(subgraph(digraph(G_filt), ...
                curr_component)));
            try
                [S, S_low, clusters, Gs] = ...
                    rate_distortion_upper_info_new(component_subgraph, setting, num_pairs);
                compressibilities_filt(k) = mean(S(end) - S);
            catch
                compressibilities_filt(k) = NaN;
            end
        end
        C_node_reordered(i, j) = sum(compressibilities_filt, 'omitnan');
    end
end

% Edge rewired networks

C_edge_rewired = zeros(num_iters, n);
DoF_edge_rewired = zeros(num_iters, n);

for i = 1:num_iters
    curr_weighted_G = edges_rewired_weighted(:, :, i);
    curr_weighted_G(curr_weighted_G == Inf) = 0;
    curr_weighted_G(curr_weighted_G ~= 0) = 1;
    curr_weighted_G(1:n+1:end) = 0;
    curr_G = curr_weighted_G;
    fprintf('Iteration %d of %d.\n', i, 10);
    for j = 1:n
        G_filt = curr_G(1:j, 1:j);
        [~, num_nodes, num_edges] = density_und(G_filt);
        DoF_edge_rewired(i, j) = (2 * num_nodes) - num_edges;
        [components, component_sizes] = conncomp(digraph(G_filt), 'Type', 'Weak');
        num_components = length(unique(components));
        compressibilities_filt = zeros(1, num_components);
        for k = 1:num_components
            curr_component = components == k;
            component_subgraph = full(adjacency(subgraph(digraph(G_filt), ...
                curr_component)));
            try
                [S, S_low, clusters, Gs] = ...
                    rate_distortion_upper_info_new(component_subgraph, setting, num_pairs);
                compressibilities_filt(k) = mean(S(end) - S);
            catch
                compressibilities_filt(k) = NaN;
            end
        end
        C_edge_rewired(i, j) = sum(compressibilities_filt, 'omitnan');
    end
end

%% Plot

figure;
hold on
plot(1:n, DoF, 'Color', [0, 0, 0], ...
    'LineWidth', 2);
% plot(1:n, mean(DoF_edge_rewired), 'LineWidth', 2, ...
%     'Color', [0.8500, 0.3250, 0.0980]);
plot(1:n, mean(DoF_node_reordered), 'LineWidth', 2, ...
    'Color', [0.9290, 0.6940, 0.1250]);
hold off
xlabel('Nodes', 'FontSize', 20);
ylabel('DoF', 'FontSize', 20);
% legend('Original', 'Edges Rewired', 'Nodes Reordered', ...
%     'Location', 'NW');
legend('Original', 'Nodes Reordered', ...
    'Location', 'NW');
title(strcat({'Subj. '}, string(subj_ID), {': DoF'}));
prettify

figure;
hold on
plot(1:n, C, 'Color', [0, 0, 0], ...
    'LineWidth', 2);
plot(1:n, mean(C_edge_rewired), 'LineWidth', 2, ...
    'Color', [0.8500, 0.3250, 0.0980]);
plot(1:n, mean(C_node_reordered), 'LineWidth', 2, ...
    'Color', [0.9290, 0.6940, 0.1250]);
hold off
xlabel('Nodes', 'FontSize', 20);
ylabel('Compressibility', 'FontSize', 20);
legend('Original', 'Edges Rewired', 'Nodes Reordered', ...
    'Location', 'NW');
title(strcat({'Subj. '}, string(subj_ID), {': Compressibility'}));
prettify

% figure;
% plot(1:n, DoF, 'Color', [0, 0, 0], ...
%     'LineWidth', 2);
% 
% figure;
% plot(1:n, mean(DoF_edge_rewired), 'LineWidth', 2, ...
%     'Color', [0.8500, 0.3250, 0.0980]);
% 
% figure;
% plot(1:n, mean(DoF_node_reordered), 'LineWidth', 2, ...
%     'Color', [0.9290, 0.6940, 0.1250]);
