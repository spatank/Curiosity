clc; close all; clear;

%% Preparatory

base_path = '/Volumes/My Passport/Curiosity/';
addpath(genpath(fullfile(base_path, 'Helper')))
addpath(genpath('/Users/sppatankar/Documents/MATLAB/humanStructureFunction'))
data_path = fullfile(base_path, 'v4/Data/KNOT/Processed_Bettis/');

subj_ID = 106;
load(strcat(data_path, 'subj_', string(subj_ID), '_bettis.mat'));

%% Preparatory

betti_0_edges_rewired_mean = mean(bettis_0_edges_rewired);
betti_1_edges_rewired_mean = mean(bettis_1_edges_rewired);
betti_2_edges_rewired_mean = mean(bettis_2_edges_rewired);

betti_0_nodes_reordered_mean =  mean(bettis_0_nodes_reordered);
betti_1_nodes_reordered_mean =  mean(bettis_1_nodes_reordered);
betti_2_nodes_reordered_mean =  mean(bettis_2_nodes_reordered);

% betti_0_edges_rewired_mean = bettis_edges_rewired_bettis_0;
% betti_1_edges_rewired_mean = bettis_edges_rewired_bettis_1;
% betti_2_edges_rewired_mean = bettis_edges_rewired_bettis_2;
% 
% betti_0_nodes_reordered_mean =  bettis_nodes_reordered_bettis_0;
% betti_1_nodes_reordered_mean =  bettis_nodes_reordered_bettis_1;
% betti_2_nodes_reordered_mean =  bettis_nodes_reordered_bettis_2;

n = size(bettis_orig, 2);

figure;
hold on
plot(1:n, bettis_orig(1, :), 'Color', [0, 0, 0], ...
    'LineWidth', 2);
plot(1:n, betti_0_edges_rewired_mean, 'LineWidth', 2, ...
    'Color', [0.8500, 0.3250, 0.0980]);
plot(1:n, betti_0_nodes_reordered_mean, 'LineWidth', 2, ...
    'Color', [0.9290, 0.6940, 0.1250]);
hold off
xlabel('Nodes', 'FontSize', 20);
ylabel('Dimension 0 Betti Number', 'FontSize', 20);
legend('Original', 'Edges Rewired', 'Nodes Reordered', ...
    'Location', 'NW');
title(strcat({'Subj. '}, string(subj_ID), {': Dimension 0 Betti Number'}));
prettify

figure;
hold on
plot(1:n, bettis_orig(2, :), 'Color', [0, 0, 0], ...
    'LineWidth', 2);
plot(1:n, betti_1_edges_rewired_mean, 'LineWidth', 2, ...
    'Color', [0.8500, 0.3250, 0.0980]);
plot(1:n, betti_1_nodes_reordered_mean, 'LineWidth', 2, ...
    'Color', [0.9290, 0.6940, 0.1250]);
hold off
xlabel('Nodes', 'FontSize', 20);
ylabel('Dimension 1 Betti Number', 'FontSize', 20);
legend('Original', 'Edges Rewired', 'Nodes Reordered', ...
    'Location', 'NW');
title(strcat({'Subj. '}, string(subj_ID), {': Dimension 1 Betti Number'}));
prettify

figure;
hold on
plot(1:n, bettis_orig(3, :), 'Color', [0, 0, 0], ...
    'LineWidth', 2);
plot(1:n, betti_2_edges_rewired_mean, 'LineWidth', 2, ...
    'Color', [0.8500, 0.3250, 0.0980]);
plot(1:n, betti_2_nodes_reordered_mean, 'LineWidth', 2, ...
    'Color', [0.9290, 0.6940, 0.1250]);
hold off
xlabel('Nodes', 'FontSize', 20);
ylabel('Dimension 2 Betti Number', 'FontSize', 20);
legend('Original', 'Edges Rewired', 'Nodes Reordered', ...
    'Location', 'NW');
title(strcat({'Subj. '}, string(subj_ID), {': Dimension 2 Betti Number'}));
prettify