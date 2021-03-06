clc; close all; clear;

%% Compressibility for Growing Knowledge Networks

base_path = '/Volumes/My Passport/Curiosity/';
addpath(genpath(fullfile(base_path, 'Helper')))
addpath(genpath(fullfile(base_path, 'Data')))
data_path = fullfile(base_path, 'v2/Data/Wiki/Wiki_processed_Eirene/');


topic = 'molecular_biology';
data_path = fullfile(data_path, strcat(topic, '.mat'));
load(data_path);

setting = 7;
num_pairs = 100;

%% Compute Compressibility

compressibility = zeros(1, n);
degrees_of_freedom = zeros(1, n);

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
    [num_nodes, num_edges] = network_size(largest_G);
    degrees_of_freedom(j) = (2 * num_nodes) - num_edges;
    try
        [S, S_low, clusters, Gs] = rate_distortion_upper_info_new(largest_G, setting, num_pairs);
        compressibility(j) = mean(S(end) - S);
    catch
        compressibility(j) = NaN;
    end
end


%% Plot Time vs. Compressibility 

% figure;
% plot(1:n, compressibility, 'LineWidth', 2, ...
%     'Color', [0, 0, 0]);
% xlabel('Node', 'FontSize', 20);
% ylabel('Compressibiility', 'FontSize', 20);
% title('Software Engineering', 'FontSize', 20);
% prettify

betti_dim_1_x = betti_curves{1, 1}(1:end - 1, 1);
betti_dim_1_y = betti_curves{1, 1}(1:end - 1, 2);

betti_dim_2_x = betti_curves{2, 1}(1:end - 1, 1);
betti_dim_2_y = betti_curves{2, 1}(1:end - 1, 2);

% figure;
% hold on
% plot(betti_dim_1_x, betti_dim_1_y, 'LineWidth', 2);
% plot(betti_dim_2_x, betti_dim_2_y, 'LineWidth', 2);
% xlabel('Node', 'FontSize', 20);
% ylabel('Betti Number', 'FontSize', 20);
% title('Optics', 'FontSize', 20);
% legend('Dimension 1', 'Dimension 2', 'Location', 'NorthWest');
% prettify

zscore_nan = @(x) bsxfun(@rdivide, bsxfun(@minus, x, mean(x, 'omitnan')), ...
    std(x, 'omitnan'));

figure;
hold on
plot(1:n, zscore_nan(compressibility), 'LineWidth', 2, ...
    'Color', [0, 0, 0]);
plot(1:n, zscore(degrees_of_freedom), 'LineWidth', 2, ...
    'Color', [0.8500, 0.3250, 0.0980]);
plot(betti_dim_1_x, zscore(betti_dim_1_y), 'LineWidth', 2, ...
    'Color', [0.9290, 0.6940, 0.1250]);
% plot(betti_dim_2_x, zscore(betti_dim_2_y), 'LineWidth', 2);
xlabel('Node', 'FontSize', 20);
ylabel('Z-Score', 'FontSize', 20);
title('Molecular Biology', 'FontSize', 20);
legend('Compressibility', '2*n - e', 'Dimension 1', 'Location', 'NorthWest');
prettify


data_1 = betti_dim_1_y;
data_2 = compressibility;
data_1_new = interp1(1:numel(data_1), data_1, linspace(1, numel(data_1), numel(data_2)));
[r, p] = corr(data_2', data_1_new', 'rows', 'complete', 'type', 'Spearman');

%% Curve Fitting

n = double(n);

x_1 = double(2:n); 
y_1 = zscore_nan(compressibility);
y_1 = y_1(2:end);
% [p1, S1] = polyfit(x_1, y_1, 2); % degree 2 polynomial
S1 = polyfitn(x_1, y_1, 2); % degree 2 polynomial
x_1_fit = linspace(1, n, 1000);
% y_1_fit = polyval(p1, x_1_fit);
y_1_fit = polyval(S1.Coefficients, x_1_fit);

x_2 = betti_dim_1_x;
y_2 = zscore(betti_dim_1_y);
% [p2, S2] = polyfit(x_2, y_2, 2); % degree 2 polynomial
S2 = polyfitn(x_2, y_2, 2); % degree 2 polynomial
x_2_fit = linspace(min(0), max(x_2), 1000);
% y_2_fit = polyval(p2, x_2_fit);
y_2_fit = polyval(S2.Coefficients, x_2_fit);

x_3 = betti_dim_2_x;
y_3 = zscore(betti_dim_2_y);
% [p3, S3] = polyfit(x_3, y_3, 2); % degree 2 polynomial
S3 = polyfitn(x_3, y_3, 2); % degree 2 polynomial
x_3_fit = linspace(min(0), max(x_3), 1000);
% y_3_fit = polyval(p3, x_3_fit);
y_3_fit = polyval(S3.Coefficients, x_3_fit);


figure;
hold on
plot(2:n, y_1, 'LineWidth', 2, ...
    'LineStyle', '--', 'Color', [0, 0, 0]);
plot(x_2, y_2, 'LineWidth', 2, ...
    'LineStyle', '--');
plot(x_3, y_3, 'LineWidth', 2, ...
    'LineStyle', '--');
xlabel('Node', 'FontSize', 20);
ylabel('Z-Score', 'FontSize', 20);
title('Molecular Biology', 'FontSize', 20);
legend('Compressibility', 'Dimension 1', 'Dimension 2', 'Location', 'NorthWest');
prettify


figure;
hold on
plot(2:n, y_1, 'LineWidth', 2, ...
    'LineStyle', '--', 'Color', [0, 0, 0]);
plot(x_2, y_2, 'LineWidth', 2, ...
    'LineStyle', '--', 'Color', [0.8500, 0.3250, 0.0980]);
plot(x_3, y_3, 'LineWidth', 2, ...
    'LineStyle', '--', 'Color', [0.9290, 0.6940, 0.1250]);
plot(x_1_fit, y_1_fit, 'LineWidth', 2, ...
    'Color', [0, 0, 0]);
plot(x_2_fit, y_2_fit, 'LineWidth', 2, ...
    'Color', [0.8500, 0.3250, 0.0980]);
plot(x_3_fit, y_3_fit, 'LineWidth', 2, ...
    'Color', [0.9290, 0.6940, 0.1250]);
xlabel('Node', 'FontSize', 20);
ylabel('Z-Score', 'FontSize', 20);
title('Molecular Biology', 'FontSize', 20);
legend('Compressibility', 'Dimension 1', 'Dimension 2', 'Location', 'NorthWest');
prettify


% Plot gradients
figure;
hold on
plot(x_1_fit, gradient(y_1_fit, x_1_fit(2) - x_1_fit(1)), 'LineWidth', 2, ...
    'Color', [0, 0, 0]);
plot(x_2_fit, gradient(y_2_fit, x_2_fit(2) - x_2_fit(1)), 'LineWidth', 2, ...
    'Color', [0.8500, 0.3250, 0.0980]);
plot(x_3_fit, gradient(y_3_fit, x_2_fit(2) - x_2_fit(1)), 'LineWidth', 2, ...
    'Color', [0.9290, 0.6940, 0.1250]);
xlabel('Node', 'FontSize', 20);
ylabel('Gradient', 'FontSize', 20);
title('Molecular Biology', 'FontSize', 20);
legend('Compressibility', 'Dimension 1', 'Dimension 2', 'Location', 'NorthWest');
prettify