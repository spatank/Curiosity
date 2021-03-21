% clc; close all; clear;

%% Preparatory

base_path = '/Volumes/My Passport/Curiosity/';
addpath(genpath(fullfile(base_path, 'Helper')))
addpath(genpath('/Users/sppatankar/Documents/MATLAB/humanStructureFunction'))

topic = 'geometry';

orig_data = load(fullfile(base_path, 'v3/Data/', ...
    strcat('all_', topic, '_processed.mat')));
null_data = load(fullfile(base_path, 'v3/Data/', ...
    strcat('all_reordered_', topic, '_processed.mat')));

%% Load and preprocess data

max_size = size(null_data.compressibilities, 2); 

compressibilities_null = real(null_data.compressibilities);
degrees_of_freedom_null = null_data.degrees_of_freedom;
betti_dim_1_null = null_data.betti_dim_1;
betti_dim_2_null = null_data.betti_dim_2;
betti_dim_3_null = null_data.betti_dim_3;

% compressibilities_null = mean(compressibilities_null, 'omitnan');
% degrees_of_freedom_null = mean(degrees_of_freedom_null, 'omitnan');
% betti_dim_1_null = mean(betti_dim_1_null, 'omitnan');
% betti_dim_2_null = mean(betti_dim_2_null, 'omitnan');
% betti_dim_3_null = mean(betti_dim_3_null, 'omitnan');
% 
% compressibilities_err = std(compressibilities_null, 'omitnan');
% degrees_of_freedom_err = std(degrees_of_freedom_null, 'omitnan');
% betti_dim_1_err = std(betti_dim_1_null, 'omitnan');
% betti_dim_2_err = std(betti_dim_2_null, 'omitnan');
% betti_dim_3_err = std(betti_dim_3_null, 'omitnan');

% compressibilities_null = mean(compressibilities_null, 'omitnan');
% degrees_of_freedom_null = mean(degrees_of_freedom_null, 'omitnan');
% betti_dim_1_null = mean(betti_dim_1_null, 'omitnan');
% betti_dim_2_null = mean(betti_dim_2_null, 'omitnan');
% betti_dim_3_null = mean(betti_dim_3_null, 'omitnan');

compressibilities_orig = real(orig_data.compressibilities);
degrees_of_freedom_orig = orig_data.degrees_of_freedom;
betti_dim_1_orig = orig_data.betti_dim_1;
betti_dim_2_orig = orig_data.betti_dim_2;
betti_dim_3_orig = orig_data.betti_dim_3;


%% Make plots

figure;
hold on
plot(1:max_size, compressibilities_null, 'Color', [0.7, 0.7, 0.7], ...
    'LineWidth', 1);
plot(1:max_size, compressibilities_orig, 'LineWidth', 2, ...
    'Color', [0, 0, 0]);
hold off
xlabel('Nodes', 'FontSize', 20);
ylabel('Compressibility', 'FontSize', 20);
title('Geometry: Compressibility');
prettify

figure;
hold on
plot(1:max_size, degrees_of_freedom_null, 'Color', [0.7, 0.7, 0.7], ...
    'LineWidth', 1);
plot(1:max_size, degrees_of_freedom_orig, 'LineWidth', 2, ...
    'Color', [0, 0, 0]);
hold off
xlabel('Nodes', 'FontSize', 20);
ylabel('2n - e', 'FontSize', 20);
title('Geometry: Degrees of Freedom');
prettify

figure;
hold on
plot(1:max_size, betti_dim_1_null, 'Color', [0.7, 0.7, 0.7], ...
    'LineWidth', 1);
plot(1:max_size, betti_dim_1_orig, 'LineWidth', 2, ...
    'Color', [0, 0, 0]);
hold off
xlabel('Nodes', 'FontSize', 20);
ylabel('Dimension 1 Betti Number', 'FontSize', 20);
title('Geometry: Dimension 1 Betti Number');
prettify