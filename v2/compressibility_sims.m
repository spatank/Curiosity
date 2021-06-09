
clc; close all; clear;

base_path = '/Volumes/My Passport/Curiosity/';
addpath(genpath(fullfile(base_path, 'Helper')))
addpath(genpath(fullfile(base_path, 'v2/Data/Simulations')))
% load('persistent_homology.mat') % constant probability model
load('persistent_homology_3.mat') % proportional probability model

%% Compute compressibility

setting = 7;
num_pairs = 100;

compressibilities = zeros(iters, n);

for i = 1:iters
    fprintf('Iteration %d of %d\n', i, iters);
    G = all_weighted_Gs(:, :, i);
    for j = 1:n
        G_filt = zeros(n, n);
        G_filt(1:j, 1:j) = G(1:j, 1:j);
        G_filt(G_filt == 2 * n) = 0; % set 0 weight edges to 0
        G_filt(G_filt > 0) = 1; % binarize
        [components, component_sizes] = conncomp(digraph(G_filt), 'Type', 'Weak');
        idx = component_sizes(components) == max(component_sizes);
        largest_G = full(adjacency(subgraph(digraph(G_filt), idx)));
        try
            [S, S_low, clusters, Gs] = rate_distortion_upper_info_new(largest_G, setting, num_pairs);
            compressibilities(i, j) = mean(S(end) - S);
        catch
            compressibilities(i, j) = NaN;
        end
    end
end

%% Plot Betti + compressibility curves

% clc; close all; clear;

% load('const_prob.mat');
% load('prop_prob.mat');
% load('osc_prob.mat');

mean_compressibility = mean(compressibilities, 'omitnan');
mean_betti_dim_1 = mean(bettiCurve(:, :, 1), 'omitnan');
mean_betti_dim_2 = mean(bettiCurve(:, :, 2), 'omitnan');
mean_betti_dim_3 = mean(bettiCurve(:, :, 3), 'omitnan');

figure;
hold on
plot(1:n, compressibilities, 'LineWidth', 2, ...
    'Color', [0.7, 0.7, 0.7]);
plot(1:n, mean_compressibility, 'LineWidth', 2, ...
    'Color', [0, 0, 0]);
hold off
xlabel('Nodes', 'FontSize', 20);
ylabel('Compressibility', 'FontSize', 20);
title('Constant Probability Model', 'FontSize', 15);
prettify

figure
hold on
plot(1:n, mean_betti_dim_1(1:n), 'LineWidth', 2, 'Color', [0, 0, 1]);
plot(1:n, mean_betti_dim_2(1:n), 'LineWidth', 2, 'Color', [0, 0.5, 0]);
plot(1:n, mean_betti_dim_3(1:n), 'LineWidth', 2, 'Color', [1, 0, 0]);
xlabel('Nodes', 'FontSize', 20);
ylabel('Betti Number', 'FontSize', 20);
title('Constant Probability Model', 'FontSize', 15);
legend('dim = 1', 'dim = 2', 'dim = 3', 'Location', 'NorthWest');
hold off
prettify

figure
hold on
data = bettiCurve(:, :, 2);
plot(1:n, data(:, 1:n), 'LineWidth', 2, 'Color', [0.7, 0.7, 0.7]);
plot(1:n, mean_betti_dim_2(1:n), 'LineWidth', 2, 'Color', [0, 0, 0]);
hold off
xlabel('Nodes', 'FontSize', 20);
ylabel('Cycles', 'FontSize', 20);
title('Constant Probability Model', 'FontSize', 15);
prettify

figure;
plot(mean_compressibility, mean_betti_dim_1(1:n), 'LineWidth', 2, ...
    'Color', [0, 0, 0]);
xlabel('Compressibility', 'FontSize', 20);
ylabel('Betti Number', 'FontSize', 20);
title('Dimension 1', 'FontSize', 15);
prettify

figure;
plot(mean_compressibility, mean_betti_dim_2(1:n), 'LineWidth', 2, ...
    'Color', [0, 0, 0]);
xlabel('Compressibility', 'FontSize', 20);
ylabel('Betti Number', 'FontSize', 20);
title('Dimension 2', 'FontSize', 15);
prettify

figure;
plot(mean_compressibility, mean_betti_dim_3(1:n), 'LineWidth', 2, ...
    'Color', [0, 0, 0]);
xlabel('Compressibility', 'FontSize', 20);
ylabel('Betti Number', 'FontSize', 20);
title('Dimension 3', 'FontSize', 15);
prettify