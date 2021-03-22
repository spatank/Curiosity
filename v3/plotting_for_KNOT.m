clc; close all; clear;

%% Preparatory

base_path = '/Volumes/My Passport/Curiosity/';
addpath(genpath(fullfile(base_path, 'Helper')))
addpath(genpath('/Users/sppatankar/Documents/MATLAB/humanStructureFunction'))

orig_data = load(fullfile(base_path, 'v3/Data/all_KNOT_processed.mat'));
null_data = load(fullfile(base_path, 'v3/Data/all_reordered_KNOT_processed.mat'));

max_size = size(null_data.compressibilities, 3); 

%% Load and preprocess data

% compressibilities_null = real(null_data.compressibilities);
% degrees_of_freedom_null = null_data.degrees_of_freedom;
% betti_dim_1_null = null_data.betti_dim_1;
% betti_dim_2_null = null_data.betti_dim_2;
% betti_dim_3_null = null_data.betti_dim_3;
% 
% compressibilities_null = squeeze(mean(compressibilities_null, 2, 'omitnan'));
% degrees_of_freedom_null = squeeze(mean(degrees_of_freedom_null, 2, 'omitnan'));
% betti_dim_1_null = squeeze(mean(betti_dim_1_null, 2, 'omitnan'));
% betti_dim_2_null = squeeze(mean(betti_dim_2_null, 2, 'omitnan'));
% betti_dim_3_null = squeeze(mean(betti_dim_3_null, 2, 'omitnan'));
% 
% compressibilities_err = std(compressibilities_null, 'omitnan');
% degrees_of_freedom_err = std(degrees_of_freedom_null, 'omitnan');
% betti_dim_1_err = std(betti_dim_1_null, 'omitnan');
% betti_dim_2_err = std(betti_dim_2_null, 'omitnan');
% betti_dim_3_err = std(betti_dim_3_null, 'omitnan');
% 
% compressibilities_null = mean(compressibilities_null, 'omitnan');
% degrees_of_freedom_null = mean(degrees_of_freedom_null, 'omitnan');
% betti_dim_1_null = mean(betti_dim_1_null, 'omitnan');
% betti_dim_2_null = mean(betti_dim_2_null, 'omitnan');
% betti_dim_3_null = mean(betti_dim_3_null, 'omitnan');
% 
% compressibilities_orig = mean(real(orig_data.compressibilities), 'omitnan');
% degrees_of_freedom_orig = mean(orig_data.degrees_of_freedom, 'omitnan');
% betti_dim_1_orig = mean(orig_data.betti_dim_1, 'omitnan');
% betti_dim_2_orig = mean(orig_data.betti_dim_2, 'omitnan');
% betti_dim_3_orig = mean(orig_data.betti_dim_3, 'omitnan');

%% Make plots for individual subjects

subj_ID = 1; % 26, 29
subj_orig_compressibilities = orig_data.compressibilities(subj_ID, :);
subj_orig_degrees_of_freedom = orig_data.degrees_of_freedom(subj_ID, :);
subj_orig_betti_dim_1 = orig_data.betti_dim_1(subj_ID, :);

% Be careful here! One subject does not have Betti curves!
subj_null_compressibilities = squeeze(null_data.compressibilities(subj_ID, :, :));
subj_null_degrees_of_freedom = squeeze(null_data.degrees_of_freedom(subj_ID, :, :));
subj_null_betti_dim_1 = squeeze(null_data.betti_dim_1(subj_ID, :, :));

figure;
hold on
plot(1:max_size, subj_null_compressibilities, 'Color', [0.7, 0.7, 0.7], ...
    'LineWidth', 1);
plot(1:max_size, subj_orig_compressibilities, 'LineWidth', 2, ...
    'Color', [0, 0, 0]);
hold off
xlabel('Nodes', 'FontSize', 20);
ylabel('Compressibility', 'FontSize', 20);
title(strcat({'Idx '}, string(subj_ID), {': Compressibility'}));
prettify

figure;
hold on
plot(1:max_size, subj_null_degrees_of_freedom, 'Color', [0.7, 0.7, 0.7], ...
    'LineWidth', 1);
plot(1:max_size, subj_orig_degrees_of_freedom, 'LineWidth', 2, ...
    'Color', [0, 0, 0]);
hold off
xlabel('Nodes', 'FontSize', 20);
ylabel('2n - e', 'FontSize', 20);
title(strcat({'Idx '}, string(subj_ID), {': Degrees of Freedom'}));
prettify

figure;
hold on
plot(1:max_size, subj_null_betti_dim_1, 'Color', [0.7, 0.7, 0.7], ...
    'LineWidth', 1);
plot(1:max_size, subj_orig_betti_dim_1, 'LineWidth', 2, ...
    'Color', [0, 0, 0]);
hold off
xlabel('Nodes', 'FontSize', 20);
ylabel('Dimension 1 Betti Number', 'FontSize', 20);
title(strcat({'Idx '}, string(subj_ID), {': Dimension 1 Betti Number'}));
prettify

%% Make averaged plots

figure;
hold on
plot(1:max_size, compressibilities_null, 'Color', [0.7, 0.7, 0.7], ...
    'LineWidth', 1);
plot(1:max_size, compressibilities_orig, 'LineWidth', 2, ...
    'Color', [0, 0, 0]);
hold off
xlabel('Nodes', 'FontSize', 20);
ylabel('Compressibility', 'FontSize', 20);
title('KNOT: Compressibility');
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
title('KNOT: Degrees of Freedom');
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
title('KNOT: Dimension 1 Betti Number');
prettify