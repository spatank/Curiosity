clc; clear;

%% Load processed data

base_path = '/Volumes/My Passport/Curiosity/';
addpath(genpath(fullfile(base_path, 'Helper')))
data_path = fullfile(base_path, 'v7/Data/Simulations/Processed/');
load(fullfile(data_path, 'C_clust_coef/simulated_nets_C_clust_coef.mat'));
load(fullfile(data_path, 'Betti/simulated_nets_bettis.mat'));
load(fullfile(data_path, 'Mech/simulated_nets_mech.mat'));

n = size(C_CP, 2);
num_iters = size(C_CP, 1);

%% Plot

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

% figure;
% hold on
% plot(1:n, C_PP, 'Color', [0.7, 0.7, 0.7], ...
%     'LineWidth', 1);
% plot(1:n, mean(C_PP, 'omitnan'), 'LineWidth', 2, ...
%     'Color', [0, 0, 0]);
% hold off
% xlabel('Nodes', 'FontSize', 20);
% ylabel('Compressibility', 'FontSize', 20);
% title('Proportional Probability Model', 'FontSize', 20);
% prettify
% 
% figure;
% hold on
% plot(1:n, clust_coef_PP, 'Color', [0.7, 0.7, 0.7], ...
%     'LineWidth', 1);
% plot(1:n, mean(clust_coef_PP, 'omitnan'), 'LineWidth', 2, ...
%     'Color', [0, 0, 0]);
% hold off
% xlabel('Nodes', 'FontSize', 20);
% ylabel('Mean Clustering Coefficient', 'FontSize', 20);
% title('Proportional Probability Model', 'FontSize', 20);
% prettify

figure;
hold on
plot(1:n, bettis_1_CP, 'Color', [0.7, 0.7, 0.7], ...
    'LineWidth', 1);
plot(1:n, mean(bettis_1_CP, 'omitnan'), 'LineWidth', 2, ...
    'Color', [0, 0, 0]);
hold off
xlabel('Nodes', 'FontSize', 20);
ylabel('Cycles', 'FontSize', 20);
title('Constant Probability Model', 'FontSize', 20);
prettify

% figure;
% hold on
% plot(1:n, bettis_1_PP, 'Color', [0.7, 0.7, 0.7], ...
%     'LineWidth', 1);
% plot(1:n, mean(bettis_1_PP, 'omitnan'), 'LineWidth', 2, ...
%     'Color', [0, 0, 0]);
% hold off
% xlabel('Nodes', 'FontSize', 20);
% ylabel('Cycles', 'FontSize', 20);
% title('Proportional Probability Model', 'FontSize', 20);
% prettify

figure;
hold on
plot(1:n, d_const, 'LineWidth', 2, ...
    'Color', [0.7, 0.7, 0.7]);
plot(1:n, mean(d_const, 'omitnan'), 'LineWidth', 2, ...
    'Color', [0, 0, 0]);
hold off
xlabel('Nodes', 'FontSize', 20);
ylabel('Embedding Dimensionality', 'FontSize', 20);
title('Constant Probability Model', 'FontSize', 20);
prettify

% figure;
% hold on
% plot(1:n, DoF_const, 'LineWidth', 2, ...
%     'Color', [0.7, 0.7, 0.7]);
% plot(1:n, mean(DoF_const, 'omitnan'), 'LineWidth', 2, ...
%     'Color', [0, 0, 0]);
% hold off
% xlabel('Nodes', 'FontSize', 20);
% ylabel('DoF', 'FontSize', 20);
% title('Constant Probability Model', 'FontSize', 20);
% prettify
% 
% figure;
% hold on
% plot(1:n, rigid_const, 'LineWidth', 2, ...
%     'Color', [0.7, 0.7, 0.7]);
% plot(1:n, mean(rigid_const, 'omitnan'), 'LineWidth', 2, ...
%     'Color', [0, 0, 0]);
% hold off
% xlabel('Nodes', 'FontSize', 20);
% ylabel('Rigid Body DoF', 'FontSize', 20);
% title('Constant Probability Model', 'FontSize', 20);
% prettify

figure;
hold on
plot(1:n, conform_const, 'LineWidth', 2, ...
    'Color', [0.7, 0.7, 0.7]);
plot(1:n, mean(conform_const, 'omitnan'), 'LineWidth', 2, ...
    'Color', [0, 0, 0]);
hold off
xlabel('Nodes', 'FontSize', 20);
ylabel('Conformable DoF', 'FontSize', 20);
title('Constant Probability Model', 'FontSize', 20);
prettify

% figure;
% hold on
% plot(1:n, d_prop, 'LineWidth', 2, ...
%     'Color', [0.7, 0.7, 0.7]);
% plot(1:n, mean(d_prop, 'omitnan'), 'LineWidth', 2, ...
%     'Color', [0, 0, 0]);
% hold off
% xlabel('Nodes', 'FontSize', 20);
% ylabel('Dimensionality', 'FontSize', 20);
% title('Proportional Probability Model', 'FontSize', 20);
% prettify
% 
% figure;
% hold on
% plot(1:n, DoF_prop, 'LineWidth', 2, ...
%     'Color', [0.7, 0.7, 0.7]);
% plot(1:n, mean(DoF_prop, 'omitnan'), 'LineWidth', 2, ...
%     'Color', [0, 0, 0]);
% hold off
% xlabel('Nodes', 'FontSize', 20);
% ylabel('DoF', 'FontSize', 20);
% title('Proportional Probability Model', 'FontSize', 20);
% prettify
% 
% figure;
% hold on
% plot(1:n, rigid_prop, 'LineWidth', 2, ...
%     'Color', [0.7, 0.7, 0.7]);
% plot(1:n, mean(rigid_prop, 'omitnan'), 'LineWidth', 2, ...
%     'Color', [0, 0, 0]);
% hold off
% xlabel('Nodes', 'FontSize', 20);
% ylabel('Rigid Body DoF', 'FontSize', 20);
% title('Proportional Probability Model', 'FontSize', 20);
% prettify
% 
% figure;
% hold on
% plot(1:n, conform_prop, 'LineWidth', 2, ...
%     'Color', [0.7, 0.7, 0.7]);
% plot(1:n, mean(conform_prop, 'omitnan'), 'LineWidth', 2, ...
%     'Color', [0, 0, 0]);
% hold off
% xlabel('Nodes', 'FontSize', 20);
% ylabel('Conformable DoF', 'FontSize', 20);
% title('Proportional Probability Model', 'FontSize', 20);
% prettify