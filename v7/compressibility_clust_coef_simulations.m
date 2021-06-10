clc; close all; clear;

%% Preparatory

base_path = '/Volumes/My Passport/Curiosity/';
addpath(genpath(fullfile(base_path, 'Helper')))
addpath(genpath('/Users/sppatankar/Documents/MATLAB/humanStructureFunction'))
load(fullfile(base_path, 'v7/Data/Simulations/Preprocessed/simulated_nets.mat'));

setting = 7;
num_pairs = 100;

n = 75; % number of nodes in simulated networks
num_iters = 50; % for network models

%% Compile metrics of interest 

% constant probability model
clust_coef_CP = NaN(num_iters, n);
C_CP = NaN(num_iters, n);

% proportional probability model
clust_coef_PP = NaN(num_iters, n);
C_PP = NaN(num_iters, n);

for iter = 1:num_iters
        
    fprintf('Iteration %d of %d.\n', iter, num_iters);
    
    G = const_prob_nets(:, :, iter);
        
    for i = 1:n
        % fprintf('Nodes %d of %d.\n', i, n);
        G_filt = G(1:i, 1:i);
        % remove nodes of degree zero
        node_degree = sum(G_filt);
        keep_nodes = find(node_degree ~= 0);
        G_filt = G_filt(keep_nodes, keep_nodes);
        if isempty(G_filt)
            continue
        end
        clust_coef_CP(iter, i) = mean(clustering_coef_bu(G_filt));
        try
            [S, S_low, clusters, Gs] = ...
                rate_distortion_upper_info(G_filt, setting, num_pairs);
            C_CP(iter, i) = mean(S(end) - S);
        catch
            fprintf('Error caught; size. %d\n', i);
            C_CP(iter, i) = NaN;
        end
    end
    
    G = prop_prob_nets(:, :, iter);
        
    for i = 1:n
        % fprintf('Nodes %d of %d.\n', i, n);
        G_filt = G(1:i, 1:i);
        % remove nodes of degree zero
        node_degree = sum(G_filt);
        keep_nodes = find(node_degree ~= 0);
        G_filt = G_filt(keep_nodes, keep_nodes);
        if isempty(G_filt)
            continue
        end
        clust_coef_PP(iter, i) = mean(clustering_coef_bu(G_filt));
        try
            [S, S_low, clusters, Gs] = ...
                rate_distortion_upper_info(G_filt, setting, num_pairs);
            C_PP(iter, i) = mean(S(end) - S);
        catch
            fprintf('Error caught; size. %d\n', i);
            C_PP(iter, i) = NaN;
        end
    end

end

% Save variables of interest    
filename = 'simulated_nets_C_clust_coef.mat';
save(filename, 'clust_coef_CP', 'C_CP', 'clust_coef_PP', 'C_PP');

%% Diagnostic plotting

close all;

figure;
hold on
plot(1:n, C_CP, 'Color', [0.7, 0.7, 0.7], ...
    'LineWidth', 1);
plot(1:n, mean(C_CP, 'omitnan'), 'LineWidth', 2, ...
    'Color', [0, 0, 0]);
hold off
xlabel('Nodes', 'FontSize', 20);
ylabel('Compressibility', 'FontSize', 20);
title('Constant Probability Model', 'FontSize', 20);
prettify

figure;
hold on
plot(1:n, clust_coef_CP, 'Color', [0.7, 0.7, 0.7], ...
    'LineWidth', 1);
plot(1:n, mean(clust_coef_CP, 'omitnan'), 'LineWidth', 2, ...
    'Color', [0, 0, 0]);
hold off
xlabel('Nodes', 'FontSize', 20);
ylabel('Mean Clustering Coefficient', 'FontSize', 20);
title('Constant Probability Model', 'FontSize', 20);
prettify

figure;
hold on
plot(1:n, C_PP, 'Color', [0.7, 0.7, 0.7], ...
    'LineWidth', 1);
plot(1:n, mean(C_PP, 'omitnan'), 'LineWidth', 2, ...
    'Color', [0, 0, 0]);
hold off
xlabel('Nodes', 'FontSize', 20);
ylabel('Compressibility', 'FontSize', 20);
title('Proportional Probability Model', 'FontSize', 20);
prettify

figure;
hold on
plot(1:n, clust_coef_PP, 'Color', [0.7, 0.7, 0.7], ...
    'LineWidth', 1);
plot(1:n, mean(clust_coef_PP, 'omitnan'), 'LineWidth', 2, ...
    'Color', [0, 0, 0]);
hold off
xlabel('Nodes', 'FontSize', 20);
ylabel('Mean Clustering Coefficient', 'FontSize', 20);
title('Proportional Probability Model', 'FontSize', 20);
prettify