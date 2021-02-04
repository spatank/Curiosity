clc; close all; clear;

%% Compressibility for Growing Knowledge Networks

base_path = '/Volumes/My Passport/Curiosity/';
addpath(genpath(fullfile(base_path, 'Helper')))
data_path = fullfile(base_path, '/Data/KNOT_processed_Eirene/');
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
compressibilities = NaN(length(files), max_size);

% for i = 1:length(files)
for i = 1:1
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
        try
            [S, S_low, clusters, Gs] = rate_distortion_upper_info_new(largest_G, setting, num_pairs);
            compressibilities(i, j) = mean(S(end) - S);
        catch
            compressibilities(i, j) = NaN;
        end
    end
end

%% Plot Time vs. Compressibility 

clc; close all;

subj_interest = 1;

nodes = 1:length(compressibilities);
compressibility = smoothdata(compressibilities(subj_interest, 1:n), ...
    'movmean', 50);
betti_dim_1 = betti_curves{1, 1}(1:n, 2);
betti_dim_2 = betti_curves{2, 1}(1:n, 2);
betti_dim_3 = betti_curves{3, 1}(1:n, 2);

figure;
plot(1:n, compressibility, 'LineWidth', 2, ...
    'Color', [0, 0, 0]);
xlabel('Nodes', 'FontSize', 20);
ylabel('Compressibility', 'FontSize', 20);
prettify

figure
hold on
plot(1:n, betti_dim_1(1:n), 'LineWidth', 2, 'Color', [0, 0, 1]);
plot(1:n, betti_dim_2(1:n), 'LineWidth', 2, 'Color', [0, 0.5, 0]);
plot(1:n, betti_dim_3(1:n), 'LineWidth', 2, 'Color', [1, 0, 0]);
xlabel('Nodes', 'FontSize', 20);
ylabel('Betti Number', 'FontSize', 20);
legend('dim = 1', 'dim = 2', 'dim = 2', 'Location', 'NorthWest');
hold off
prettify

figure;
plot(compressibility, betti_dim_1(1:n), 'LineWidth', 2, ...
    'Color', [0, 0, 0]);
xlabel('Compressibility', 'FontSize', 20);
ylabel('Betti Number', 'FontSize', 20);
title('Dimension 1', 'FontSize', 15);
prettify

figure;
plot(compressibility, betti_dim_2(1:n), 'LineWidth', 2, ...
    'Color', [0, 0, 0]);
xlabel('Compressibility', 'FontSize', 20);
ylabel('Betti Number', 'FontSize', 20);
title('Dimension 2', 'FontSize', 15);
prettify

figure;
plot(compressibility, betti_dim_3(1:n), 'LineWidth', 2, ...
    'Color', [0, 0, 0]);
xlabel('Compressibility', 'FontSize', 20);
ylabel('Betti Number', 'FontSize', 20);
title('Dimension 3', 'FontSize', 15);
prettify
