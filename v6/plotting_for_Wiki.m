clc; close all; clear;

%% Load processed data

base_path = '/Volumes/My Passport/Curiosity/';
addpath(genpath(fullfile(base_path, 'Helper')))
data_path = '/Volumes/My Passport/Curiosity/v6/Data/Wiki/';
topic = 'molecular_biology';
load(strcat(data_path, 'Processed/C_DoF/', topic, '_C_DoF.mat'))
load(strcat(data_path, 'Processed/Betti/', topic, '_bettis.mat'))
load(strcat(data_path, 'Processed/d/', topic, '_d.mat'))
num_iters = size(C_edge_rewired, 1);

%% Plot

figure;
hold on
plot(1:n, DoF, 'Color', [0, 0, 0], ...
    'LineWidth', 2);
errorbar(1:n, mean(DoF_edge_rewired, 'omitnan'), std(DoF_edge_rewired, 'omitnan')/num_iters, ...
    'MarkerSize', 2, 'LineWidth', 1, 'Color', [0.8500, 0.3250, 0.0980], ...
    'MarkerEdgeColor',[0.8500, 0.3250, 0.0980], ...
    'MarkerFaceColor',[0.8500, 0.3250, 0.0980]);
errorbar(1:n, mean(DoF_latticized, 'omitnan'), std(DoF_latticized, 'omitnan')/num_iters, ...
    'MarkerSize', 2, 'LineWidth', 1, 'Color', [0.9290, 0.6940, 0.1250], ...
    'MarkerEdgeColor',[0.9290, 0.6940, 0.1250], ...
    'MarkerFaceColor',[0.9290, 0.6940, 0.1250]);
hold off
xlabel('Nodes', 'FontSize', 20);
ylabel('DoF', 'FontSize', 20);
legend('Original', 'Edges Rewired', 'Latticized', ...
    'Location', 'NW');
title('Molecular Biology');
prettify

figure;
hold on
plot(1:n, C, 'Color', [0, 0, 0], ...
    'LineWidth', 2);
errorbar(1:n, mean(C_edge_rewired, 'omitnan'), std(C_edge_rewired, 'omitnan')/num_iters, ...
    'MarkerSize', 2, 'LineWidth', 1, 'Color', [0.8500, 0.3250, 0.0980], ...
    'MarkerEdgeColor',[0.8500, 0.3250, 0.0980], ...
    'MarkerFaceColor',[0.8500, 0.3250, 0.0980]);
errorbar(1:n, mean(C_latticized, 'omitnan'), std(C_latticized, 'omitnan')/num_iters, ...
    'MarkerSize', 2, 'LineWidth', 1, 'Color', [0.9290, 0.6940, 0.1250], ...
    'MarkerEdgeColor',[0.9290, 0.6940, 0.1250], ...
    'MarkerFaceColor',[0.9290, 0.6940, 0.1250]);
hold off
xlabel('Nodes', 'FontSize', 20);
ylabel('Compressibility', 'FontSize', 20);
legend('Original', 'Edges Rewired', 'Latticized', ...
    'Location', 'NW');
title('Molecular Biology');
prettify

figure;
hold on
plot(1:n, bettis_orig(2, :), 'Color', [0, 0, 0], ...
    'LineWidth', 2);
errorbar(1:n, mean(bettis_1_edges_rewired, 'omitnan'), std(bettis_1_edges_rewired, 'omitnan')/num_iters, ...
    'MarkerSize', 2, 'LineWidth', 1, 'Color', [0.8500, 0.3250, 0.0980], ...
    'MarkerEdgeColor',[0.8500, 0.3250, 0.0980], ...
    'MarkerFaceColor',[0.8500, 0.3250, 0.0980]);
errorbar(1:n, mean(bettis_1_latticized, 'omitnan'), std(bettis_1_latticized, 'omitnan')/num_iters, ...
    'MarkerSize', 2, 'LineWidth', 1, 'Color', [0.9290, 0.6940, 0.1250], ...
    'MarkerEdgeColor',[0.9290, 0.6940, 0.1250], ...
    'MarkerFaceColor',[0.9290, 0.6940, 0.1250]);
hold off
xlabel('Nodes', 'FontSize', 20);
ylabel('Cycles', 'FontSize', 20);
legend('Original', 'Edges Rewired', 'Latticized', ...
    'Location', 'NW');
title('Molecular Biology');
prettify

figure;
hold on
plot(1:n, d, 'LineWidth', 2, ...
    'Color', [0, 0, 0]);
plot(1:n, ceil(mean(d_edge_rewired, 'omitnan')), 'LineWidth', 2, ...
    'Color', [0.8500, 0.3250, 0.0980]);
plot(1:n, ceil(mean(d_latticized, 'omitnan')), 'LineWidth', 2, ...
    'Color', [0.9290, 0.6940, 0.1250]);
hold off
xlabel('Nodes', 'FontSize', 20);
ylabel('Dimensionality', 'FontSize', 20);
legend('Original', 'Edges Rewired', 'Latticized', ...
    'Location', 'NW');
title('Molecular Biology');
prettify

zscore_nan = @(x) bsxfun(@rdivide, bsxfun(@minus, x, mean(x, 'omitnan')), ...
    std(x, 'omitnan'));

figure;
hold on
plot(1:n, zscore_nan(C), 'Color', [0, 0, 0], ...
    'LineWidth', 2);
plot(1:n, zscore_nan(bettis_orig(2, :)), 'Color', [0, 0, 0], ...
    'LineWidth', 2, 'LineStyle', '--');
hold off
xlabel('Nodes', 'FontSize', 20);
ylabel('Z-Score', 'FontSize', 20);
legend('Compressibility', 'Cycles', 'Location', 'NW');
title('Molecular Biology');
prettify