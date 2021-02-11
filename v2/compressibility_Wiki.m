clc; close all; clear;

%% Compressibility for Growing Knowledge Networks

base_path = '/Volumes/My Passport/Curiosity/';
addpath(genpath(fullfile(base_path, 'Helper')))
addpath(genpath(fullfile(base_path, 'Data')))
data_path = fullfile(base_path, 'v2/Data/Wiki/Wiki_processed_Eirene/');


topic = 'software_engineering';
data_path = fullfile(data_path, strcat(topic, '.mat'));
load(data_path);

setting = 7;
num_pairs = 100;

%% Compute Compressibility

compressibility = zeros(1, n);

G = weighted_adj;
for j = 1:n
    fprintf('%d nodes out of %d\n', j, n);
    G_filt = zeros(n, n);
    G_filt(1:j, 1:j) = G(1:j, 1:j);
    G_filt(G_filt == 2 * n) = 0; % set 0 weight edges to 0
    G_filt(G_filt > 0) = 1; % binarize
    [components, component_sizes] = conncomp(digraph(G_filt), 'Type', 'Weak');
    idx = component_sizes(components) == max(component_sizes);
    largest_G = full(adjacency(subgraph(digraph(G_filt), idx)));
    try
        [S, S_low, clusters, Gs] = rate_distortion_upper_info_new(largest_G, setting, num_pairs);
        compressibility(j) = mean(S(end) - S);
    catch
        compressibility(j) = NaN;
    end
end


%% Plot Time vs. Compressibility 

close all

figure;
plot(1:n, compressibility, 'LineWidth', 2, ...
    'Color', [0, 0, 0]);
xlabel('Node', 'FontSize', 20);
ylabel('Compressibiility', 'FontSize', 20);
title('Software Engineering', 'FontSize', 20);
prettify

betti_dim_1_x = betti_curves{1, 1}(1:end - 1, 1);
betti_dim_1_y = betti_curves{1, 1}(1:end - 1, 2);

betti_dim_2_x = betti_curves{2, 1}(1:end - 1, 1);
betti_dim_2_y = betti_curves{2, 1}(1:end - 1, 2);

figure;
hold on
plot(betti_dim_1_x, betti_dim_1_y, 'LineWidth', 2);
plot(betti_dim_2_x, betti_dim_2_y, 'LineWidth', 2);
xlabel('Node', 'FontSize', 20);
ylabel('Betti Number', 'FontSize', 20);
title('Software Engineering', 'FontSize', 20);
legend('Dimension 1', 'Dimension 2', 'Location', 'NorthWest');
prettify

zscore_nan = @(x) bsxfun(@rdivide, bsxfun(@minus, x, mean(x, 'omitnan')), ...
    std(x, 'omitnan'));

figure;
hold on
plot(1:n, zscore_nan(compressibility), 'LineWidth', 2, ...
    'Color', [0, 0, 0]);
plot(betti_dim_1_x, zscore(betti_dim_1_y), 'LineWidth', 2);
plot(betti_dim_2_x, zscore(betti_dim_2_y), 'LineWidth', 2);
xlabel('Node', 'FontSize', 20);
ylabel('Z-Score', 'FontSize', 20);
title('Software Engineering', 'FontSize', 20);
legend('Compressibility', 'Dimension 1', 'Dimension 2', 'Location', 'NorthWest');
prettify
    

data_1 = betti_dim_1_y;
data_2 = compressibility;
data_1_new = interp1(1:numel(data_1), data_1, linspace(1, numel(data_1), numel(data_2)));
[r, p] = corr(data_2', data_1_new', 'rows', 'complete', 'type', 'Spearman');