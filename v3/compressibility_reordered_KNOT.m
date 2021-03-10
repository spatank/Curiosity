clc; close all; clear;

%% Preparatory

base_path = '/Volumes/My Passport/Curiosity/';
addpath(genpath(fullfile(base_path, 'Helper')))
addpath(genpath('/Users/sppatankar/Documents/MATLAB/humanStructureFunction'))
data_path = fullfile(base_path, 'v3/Data/reordered_KNOT_preprocessed_Eirene');
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

%% Compile metrics of interest for growing knowledge networks
 
n = length(adj);
iters = size(reordered_adjs, 3);

% initialize empty matrix of compressibility values
compressibilities = NaN(length(files), iters, max_size);
degrees_of_freedom = NaN(length(files), iters, max_size);
% betti_dim_1 = NaN(length(files), max_size, iters);
% betti_dim_2 = NaN(length(files), max_size, iters);
% betti_dim_3 = NaN(length(files), max_size, iters);

for i = 1:length(files)
    load(fullfile(data_path, files(i).name));
    % adj, reordered_adjs, reordered_nodes, reordered_weighted_Gs available
    for j = 1:iters
        fprintf('Subject %d of %d. Iteration %d or %d.\n', ...
            i, length(files), j, iters)
        G = reordered_weighted_Gs(:, :, j);
        for k = 1:n
            G_filt = zeros(n, n);
            G_filt(1:k, 1:k) = G(1:k, 1:k);
            G_filt(G_filt == 2 * n) = 0; % set 0 weight edges to 0
            G_filt(G_filt > 0) = 1; % binarize
            [components, component_sizes] = conncomp(digraph(G_filt), 'Type', 'Weak');
            idx = component_sizes(components) == max(component_sizes);
            largest_G = full(adjacency(subgraph(digraph(G_filt), idx)));
            [~, num_nodes, num_edges] = density_und(largest_G);
            degrees_of_freedom(i, j, k) = (2 * num_nodes) - num_edges;
            try
                [S, S_low, clusters, Gs] = rate_distortion_upper_info_new(largest_G, setting, num_pairs);
                compressibilities(i, j) = mean(S(end) - S);
            catch
                compressibilities(i, j) = NaN;
            end
        end
    end
end

save_string = fullfile(base_path, 'v3/Data/partial_reordered_KNOT_processed.mat');
save(save_string, 'compressibilities', 'degrees_of_freedom');

%% 

load(fullfile(base_path, 'v3/Data/all_KNOT_processed.mat')) ;

max_size = size(compressibilities, 2); 

compressibilities_raw = mean(compressibilities, 'omitnan');
compressibilities = smoothdata(compressibilities_raw, 'movmean', 25);

degrees_of_freedom_raw = mean(degrees_of_freedom, 'omitnan');
degrees_of_freedom = smoothdata(degrees_of_freedom_raw, 'movmean', 25);

betti_dim_1_raw = mean(betti_dim_1, 'omitnan');
betti_dim_1_smooth = smoothdata(betti_dim_1_raw, 'movmean', 25);
betti_dim_2_raw = mean(betti_dim_2, 'omitnan');
betti_dim_2_smooth = smoothdata(betti_dim_2_raw, 'movmean', 25);
betti_dim_3_raw = mean(betti_dim_3, 'omitnan');
betti_dim_3_smooth = smoothdata(betti_dim_3_raw, 'movmean', 25);

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

figure;
hold on
plot(1:max_size, betti_dim_1_raw, 'Color', [0.7, 0.7, 0.7], ...
    'LineWidth', 2);
plot(1:max_size, betti_dim_1_smooth, 'LineWidth', 2, ...
    'Color', [0, 0, 0]);
hold off
xlabel('Nodes', 'FontSize', 20);
ylabel('Dimension 1 Betti Number', 'FontSize', 20);
title('KNOT');
prettify

zscore_nan = @(x) bsxfun(@rdivide, bsxfun(@minus, x, mean(x, 'omitnan')), ...
    std(x, 'omitnan'));

figure;
hold on
plot(1:length(compressibilities), zscore_nan(compressibilities), 'LineWidth', 2, ...
    'Color', [0, 0, 0]);
plot(1:length(degrees_of_freedom), zscore(degrees_of_freedom), 'LineWidth', 2, ...
    'Color', [0.8500, 0.3250, 0.0980]);
plot(1:length(betti_dim_1_smooth), zscore_nan(betti_dim_1_smooth), 'LineWidth', 2, ...
    'Color', [0.9290, 0.6940, 0.1250]);
xlabel('Node', 'FontSize', 20);
ylabel('Z-Score', 'FontSize', 20);
title('KNOT', 'FontSize', 20);
legend('Compressibility', '2n - e', 'Betti Number (Dim. 1)', ...
    'Location', 'NorthWest');
prettify