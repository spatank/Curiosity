clc; close all; clear;

%% Compressibility for Growing Knowledge Networks

base_path = '/Volumes/My Passport/Curiosity/';
addpath(genpath(fullfile(base_path, 'Helper')))
data_path = fullfile(base_path, 'v3/Data/KNOT_processed_Eirene/');
files = dir(fullfile(data_path, 'subj_*.mat'));

setting = 7;
num_pairs = 100;

% how large is the largest network?
max_size = 0;
for i = 1:length(files)
    load(fullfile(data_path, files(i).name));
    if size(adj, 1) > max_size
        max_size = size(adj, 1);
    end
end

%% Compile Compressibility Values

% initialize empty matrix of compressibility values
% compressibilities = NaN(length(files), max_size);
degrees_of_freedom = NaN(length(files), max_size);
% betti_dim_1 = NaN(length(files), max_size);
% betti_dim_2 = NaN(length(files), max_size);
% betti_dim_3 = NaN(length(files), max_size);

for i = 1:length(files)
    fprintf('Subject %d of %d.\n', i, length(files))
    load(fullfile(data_path, files(i).name));
    G = weighted_adj;
    for j = 1:n
        G_filt = zeros(n, n);
        G_filt(1:j, 1:j) = G(1:j, 1:j);
        G_filt(G_filt == 2 * n) = 0; % set 0 weight edges to 0
        G_filt(G_filt > 0) = 1; % binarize
        [components, component_sizes] = conncomp(digraph(G_filt), 'Type', 'Weak');
        idx = component_sizes(components) == max(component_sizes);
        largest_G = full(adjacency(subgraph(digraph(G_filt), idx)));
        [num_nodes, num_edges] = network_size(largest_G);
        degrees_of_freedom(i, j) = (2 * num_nodes) - num_edges;
%         try
%             [S, S_low, clusters, Gs] = rate_distortion_upper_info_new(largest_G, setting, num_pairs);
%             compressibilities(i, j) = mean(S(end) - S);
%         catch
%             compressibilities(i, j) = NaN;
%         end
    end
%     betti_dim_1(i, 1:n) = betti_curves{1, 1}(1:n-1, 2)';
%     % transpose is needed for alignment reasons
%     betti_dim_2(i, 1:n) = betti_curves{2, 1}(1:n-1, 2)';
%     betti_dim_3(i, 1:n) = betti_curves{3, 1}(1:n-1, 2)';
end

% save_string = fullfile(base_path, 'v3/Data/KNOT/all_KNOT_processed.mat');
% save(save_string, 'compressibilities', 'betti_dim_1', ...
%     'betti_dim_2', 'betti_dim_3');

load('test.mat')

max_size = size(compressibilities, 2); 

nodes = 1:length(compressibilities);
% compressibilities_raw = mean(compressibilities, 'omitnan');
% compressibilities = smoothdata(compressibilities_raw, 'movmean', 25);

degrees_of_freedom_raw = mean(degrees_of_freedom, 'omitnan');
degrees_of_freedom = smoothdata(degrees_of_freedom_raw, 'movmean', 25);

figure;
hold on
plot(1:max_size, compressibilities_raw, 'Color', [0.7, 0.7, 0.7], ...
    'LineWidth', 2);
plot(1:max_size, compressibilities, 'LineWidth', 2, ...
    'Color', [0, 0, 0]);
hold off
xlabel('Nodes', 'FontSize', 20);
ylabel('Compressibility', 'FontSize', 20);
title('KNOT');
prettify

figure;
hold on
plot(1:max_size, degrees_of_freedom_raw, 'Color', [0.7, 0.7, 0.7], ...
    'LineWidth', 2);
plot(1:max_size, degrees_of_freedom, 'LineWidth', 2, ...
    'Color', [0, 0, 0]);
hold off
xlabel('Nodes', 'FontSize', 20);
ylabel('2n - e', 'FontSize', 20);
title('KNOT');
prettify

%% Plot Time vs. Compressibility 

% clc; close all; clear;
% 
% addpath(genpath('/Users/sppatankar/Documents/MATLAB/humanStructureFunction/Analysis'));
% 
% load('/Volumes/My Passport/Curiosity/v3/Data/KNOT/all_KNOT_processed.mat')
% 
% max_size = size(compressibilities, 2); 
% 
% nodes = 1:length(compressibilities);
% compressibilities_raw = mean(compressibilities, 'omitnan');
% compressibilities = smoothdata(compressibilities_raw, 'movmean', 25);
% 
% betti_dim_1_raw = mean(betti_dim_1, 'omitnan');
% betti_dim_1_smooth = smoothdata(betti_dim_1_raw, 'movmean', 25);
% betti_dim_2_raw = mean(betti_dim_2, 'omitnan');
% betti_dim_2_smooth = smoothdata(betti_dim_2_raw, 'movmean', 25);
% betti_dim_3_raw = mean(betti_dim_3, 'omitnan');
% betti_dim_3_smooth = smoothdata(betti_dim_3_raw, 'movmean', 25);


% % FIGURE OF INTEREST
% figure;
% plot(compressibilities, betti_dim_1_smooth, 'LineWidth', 2, ...
%     'Color', [0, 0, 0]);
% xlabel('Compressibility', 'FontSize', 20);
% ylabel('Betti Number (Dimension 1)', 'FontSize', 20);
% title('KNOT', 'FontSize', 15);
% prettify
% 
% [r_comp_betti_dim_1, p_comp_betti_dim_1, residuals] = ... 
%     partialcorr_with_resids(compressibilities', betti_dim_1_smooth', (1:max_size)', ...
%     'type', 'Spearman', 'rows', 'complete');
% fprintf('r = %f, p = %f. comp. ~ Betti 1, controlling for size\n', r_comp_betti_dim_1, p_comp_betti_dim_1);
% comp_resid = residuals(:, 1);
% betti_dim_1_resid = residuals(:, 2);
% coeffs = polyfit(comp_resid, betti_dim_1_resid, 1);
% x = linspace(min(comp_resid), max(comp_resid), 1000);
% y = polyval(coeffs, x);
% f = figure('color', 'w');
% hold on
% plot(comp_resid, betti_dim_1_resid,  'Color', [0, 0, 0], ...
%     'LineWidth', 2);
% plot(x, y, 'LineWidth', 2);
% xlabel('Compressibility Residual', 'FontSize', 15);
% ylabel('Betti Number Residual', 'FontSize', 15);
% title('Dimension 1 (Network Size Controlled)', ...
%     'FontSize', 15);
% prettify
% hold off

% figure;
% hold on
% plot(1:max_size, compressibilities_raw, 'Color', [0.7, 0.7, 0.7], ...
%     'LineWidth', 2);
% plot(1:max_size, compressibilities, 'LineWidth', 2, ...
%     'Color', [0, 0, 0]);
% hold off
% xlabel('Nodes', 'FontSize', 20);
% ylabel('Compressibility', 'FontSize', 20);
% title('KNOT');
% prettify

% figure;
% hold on
% plot(1:max_size, betti_dim_1_raw, 'Color', [0.7, 0.7, 0.7], ...
%     'LineWidth', 2);
% plot(1:max_size, betti_dim_1_smooth, 'LineWidth', 2, ...
%     'Color', [0, 0, 0]);
% hold off
% xlabel('Nodes', 'FontSize', 20);
% ylabel('Betti Number (Dimension 1)', 'FontSize', 20);
% title('KNOT');
% prettify
