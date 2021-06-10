clc; close all; clear;

%% Load processed data

base_path = '/Volumes/My Passport/Curiosity/';
addpath(genpath(fullfile(base_path, 'Helper')))
data_path = '/Volumes/My Passport/Curiosity/v7/Data/Wiki/';
topic = 'molecular_biology';
load(strcat(data_path, 'Processed/C_clust_coef/', topic, '_C_clust_coef.mat'))
load(strcat(data_path, 'Processed/Betti/', topic, '_bettis.mat'))
load(strcat(data_path, 'Processed/Mech/', topic, '_mech.mat'))
num_iters = size(C_edge_rewired, 1);

%% Plot

figure;
hold on
plot(1:n, C, 'Color', [0, 0, 0], ...
    'LineWidth', 1);
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
title('Molecular Biology: Compressibility', 'FontSize', 20);
prettify

figure;
hold on
plot(1:n, clust_coef, 'Color', [0, 0, 0], ...
    'LineWidth', 1);
errorbar(1:n, mean(clust_coef_edge_rewired, 'omitnan'), std(clust_coef_edge_rewired, 'omitnan')/num_iters, ...
    'MarkerSize', 2, 'LineWidth', 1, 'Color', [0.8500, 0.3250, 0.0980], ...
    'MarkerEdgeColor',[0.8500, 0.3250, 0.0980], ...
    'MarkerFaceColor',[0.8500, 0.3250, 0.0980]);
errorbar(1:n, mean(clust_coef_latticized, 'omitnan'), std(clust_coef_latticized, 'omitnan')/num_iters, ...
    'MarkerSize', 2, 'LineWidth', 1, 'Color', [0.9290, 0.6940, 0.1250], ...
    'MarkerEdgeColor',[0.9290, 0.6940, 0.1250], ...
    'MarkerFaceColor',[0.9290, 0.6940, 0.1250]);
hold off
xlabel('Nodes', 'FontSize', 20);
ylabel('Mean Clustering Coefficient', 'FontSize', 20);
legend('Original', 'Edges Rewired', 'Latticized', ...
    'Location', 'NW');
title('Molecular Biology: Clustering Coefficient', 'FontSize', 20);
prettify

figure;
hold on
plot(1:n, bettis_orig(2, :), 'Color', [0, 0, 0], ...
    'LineWidth', 1);
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
ylabel('Betti Number', 'FontSize', 20);
legend('Original', 'Edges Rewired', 'Latticized', ...
    'Location', 'NW');
title('Molecular Biology: Cycles', 'FontSize', 20);
prettify

figure;
hold on
plot(1:n, d, 'Color', [0, 0, 0], ...
    'LineWidth', 1);
errorbar(1:n, mean(d_edge_rewired, 'omitnan'), std(d_edge_rewired, 'omitnan')/num_iters, ...
    'MarkerSize', 2, 'LineWidth', 1, 'Color', [0.8500, 0.3250, 0.0980], ...
    'MarkerEdgeColor',[0.8500, 0.3250, 0.0980], ...
    'MarkerFaceColor',[0.8500, 0.3250, 0.0980]);
errorbar(1:n, mean(d_latticized, 'omitnan'), std(d_latticized, 'omitnan')/num_iters, ...
    'MarkerSize', 2, 'LineWidth', 1, 'Color', [0.9290, 0.6940, 0.1250], ...
    'MarkerEdgeColor',[0.9290, 0.6940, 0.1250], ...
    'MarkerFaceColor',[0.9290, 0.6940, 0.1250]);
hold off
xlabel('Nodes', 'FontSize', 20);
ylabel('Embedding Dimensionality', 'FontSize', 20);
legend('Original', 'Edges Rewired', 'Latticized', ...
    'Location', 'NW');
title('Molecular Biology: Embedding Dimensionality', 'FontSize', 20);
prettify

figure;
hold on
plot(1:n, DoF, 'Color', [0, 0, 0], ...
    'LineWidth', 1);
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
title('Molecular Biology: DoF', 'FontSize', 20);
prettify

figure;
hold on
plot(1:n, rigid, 'Color', [0, 0, 0], ...
    'LineWidth', 1);
errorbar(1:n, mean(rigid_edge_rewired, 'omitnan'), std(rigid_edge_rewired, 'omitnan')/num_iters, ...
    'MarkerSize', 2, 'LineWidth', 1, 'Color', [0.8500, 0.3250, 0.0980], ...
    'MarkerEdgeColor',[0.8500, 0.3250, 0.0980], ...
    'MarkerFaceColor',[0.8500, 0.3250, 0.0980]);
errorbar(1:n, mean(rigid_latticized, 'omitnan'), std(rigid_latticized, 'omitnan')/num_iters, ...
    'MarkerSize', 2, 'LineWidth', 1, 'Color', [0.9290, 0.6940, 0.1250], ...
    'MarkerEdgeColor',[0.9290, 0.6940, 0.1250], ...
    'MarkerFaceColor',[0.9290, 0.6940, 0.1250]);
hold off
xlabel('Nodes', 'FontSize', 20);
ylabel('Rigid Body DoF', 'FontSize', 20);
legend('Original', 'Edges Rewired', 'Latticized', ...
    'Location', 'NW');
title('Molecular Biology: Rigid Body DoF', 'FontSize', 20);
prettify

figure;
hold on
plot(1:n, conform, 'Color', [0, 0, 0], ...
    'LineWidth', 1);
errorbar(1:n, mean(conform_edge_rewired, 'omitnan'), std(conform_edge_rewired, 'omitnan')/num_iters, ...
    'MarkerSize', 2, 'LineWidth', 1, 'Color', [0.8500, 0.3250, 0.0980], ...
    'MarkerEdgeColor',[0.8500, 0.3250, 0.0980], ...
    'MarkerFaceColor',[0.8500, 0.3250, 0.0980]);
errorbar(1:n, mean(conform_latticized, 'omitnan'), std(conform_latticized, 'omitnan')/num_iters, ...
    'MarkerSize', 2, 'LineWidth', 1, 'Color', [0.9290, 0.6940, 0.1250], ...
    'MarkerEdgeColor',[0.9290, 0.6940, 0.1250], ...
    'MarkerFaceColor',[0.9290, 0.6940, 0.1250]);
hold off
xlabel('Nodes', 'FontSize', 20);
ylabel('Conformable DoF', 'FontSize', 20);
legend('Original', 'Edges Rewired', 'Latticized', ...
    'Location', 'NW');
title('Molecular Biology: Conformable DoF', 'FontSize', 20);
prettify