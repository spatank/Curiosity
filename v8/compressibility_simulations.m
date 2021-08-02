clc; close all; clear;

tic

%% Preparatory

base_path = '/Volumes/My Passport/Curiosity/';
addpath(genpath(fullfile(base_path, 'Helper')))
addpath(genpath('/Users/sppatankar/Documents/MATLAB/humanStructureFunction'))
load(fullfile(base_path, 'v8/Data/Simulations/simulated_nets.mat'));

setting = 7;
num_pairs = 100;

n = size(const_prob_nets, 1); % number of nodes in simulated networks
num_iters = size(const_prob_nets, 3); % for network models

%% Compile metrics of interest 

% constant probability model
clust_coef_CP = NaN(num_iters, n);
C_CP = NaN(num_iters, n);
C_norm_CP = NaN(num_iters, n);
entropy_CP = NaN(num_iters, n);


% proportional probability model
clust_coef_PP = NaN(num_iters, n);
C_PP = NaN(num_iters, n);
C_norm_PP = NaN(num_iters, n);
entropy_PP = NaN(num_iters, n);

% preferential attachment model
clust_coef_PA = NaN(num_iters, n);
C_PA = NaN(num_iters, n);
C_norm_PA = NaN(num_iters, n);
entropy_PA = NaN(num_iters, n);


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
            C_norm_CP(iter, i) = C_CP(iter, i)/S(end);
            entropy_CP(iter, i) = S(end);
        catch
            fprintf('Error caught; size. %d\n', i);
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
            C_norm_PP(iter, i) = C_PP(iter, i)/S(end);
            entropy_PP(iter, i) = S(end);
        catch
            fprintf('Error caught; size. %d\n', i);
        end
    end
    
    G = PA_nets(:, :, iter);
        
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
        clust_coef_PA(iter, i) = mean(clustering_coef_bu(G_filt));
        try
            [S, S_low, clusters, Gs] = ...
                rate_distortion_upper_info(G_filt, setting, num_pairs);
            C_PA(iter, i) = mean(S(end) - S);
            C_norm_PA(iter, i) = C_PA(iter, i)/S(end);
            entropy_PA(iter, i) = S(end);
        catch
            fprintf('Error caught; size. %d\n', i);
        end
    end

end

% Save variables of interest    
filename = 'simulated_nets_C.mat';
save(filename, 'clust_coef_CP', 'C_CP', 'C_norm_CP', 'entropy_CP', ...
    'clust_coef_PP', 'C_PP', 'C_norm_PP', 'entropy_PP', ...
    'clust_coef_PA', 'C_PA', 'C_norm_PA', 'entropy_PA');



%% Diagnostic plotting

n = 70;

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
plot(1:n, C_norm_CP, 'Color', [0.7, 0.7, 0.7], ...
    'LineWidth', 1);
plot(1:n, mean(C_norm_CP, 'omitnan'), 'LineWidth', 2, ...
    'Color', [0, 0, 0]);
hold off
xlabel('Nodes', 'FontSize', 20);
ylabel('Normalized Compressibility', 'FontSize', 20);
title('Constant Probability Model', 'FontSize', 20);
prettify

figure;
hold on
plot(1:n, entropy_CP, 'Color', [0.7, 0.7, 0.7], ...
    'LineWidth', 1);
plot(1:n, mean(entropy_CP, 'omitnan'), 'LineWidth', 2, ...
    'Color', [0, 0, 0]);
hold off
xlabel('Nodes', 'FontSize', 20);
ylabel('Entropy', 'FontSize', 20);
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
plot(1:n, C_norm_PP, 'Color', [0.7, 0.7, 0.7], ...
    'LineWidth', 1);
plot(1:n, mean(C_norm_PP, 'omitnan'), 'LineWidth', 2, ...
    'Color', [0, 0, 0]);
hold off
xlabel('Nodes', 'FontSize', 20);
ylabel('Normalized Compressibility', 'FontSize', 20);
title('Proportional Probability Model', 'FontSize', 20);
prettify

figure;
hold on
plot(1:n, entropy_PP, 'Color', [0.7, 0.7, 0.7], ...
    'LineWidth', 1);
plot(1:n, mean(entropy_PP, 'omitnan'), 'LineWidth', 2, ...
    'Color', [0, 0, 0]);
hold off
xlabel('Nodes', 'FontSize', 20);
ylabel('Entropy', 'FontSize', 20);
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

figure;
hold on
plot(1:n, C_PA, 'Color', [0.7, 0.7, 0.7], ...
    'LineWidth', 1);
plot(1:n, mean(C_PA, 'omitnan'), 'LineWidth', 2, ...
    'Color', [0, 0, 0]);
hold off
xlabel('Nodes', 'FontSize', 20);
ylabel('Compressibility', 'FontSize', 20);
title('Preferential Attachment Model', 'FontSize', 20);
prettify

figure;
hold on
plot(1:n, C_norm_PA, 'Color', [0.7, 0.7, 0.7], ...
    'LineWidth', 1);
plot(1:n, mean(C_norm_PA, 'omitnan'), 'LineWidth', 2, ...
    'Color', [0, 0, 0]);
hold off
xlabel('Nodes', 'FontSize', 20);
ylabel('Normalized Compressibility', 'FontSize', 20);
title('Preferential Attachment Model', 'FontSize', 20);
prettify

figure;
hold on
plot(1:n, entropy_PA, 'Color', [0.7, 0.7, 0.7], ...
    'LineWidth', 1);
plot(1:n, mean(entropy_PA, 'omitnan'), 'LineWidth', 2, ...
    'Color', [0, 0, 0]);
hold off
xlabel('Nodes', 'FontSize', 20);
ylabel('Entropy', 'FontSize', 20);
title('Preferential Attachment Model', 'FontSize', 20);
prettify

figure;
hold on
plot(1:n, clust_coef_PA, 'Color', [0.7, 0.7, 0.7], ...
    'LineWidth', 1);
plot(1:n, mean(clust_coef_PA, 'omitnan'), 'LineWidth', 2, ...
    'Color', [0, 0, 0]);
hold off
xlabel('Nodes', 'FontSize', 20);
ylabel('Mean Clustering Coefficient', 'FontSize', 20);
title('Preferential Attachment Model', 'FontSize', 20);
prettify