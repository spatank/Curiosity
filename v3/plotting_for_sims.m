clc; close all; clear;

%% Preparatory

base_path = '/Volumes/My Passport/Curiosity/';
addpath(genpath(fullfile(base_path, 'Helper')))
addpath(genpath('/Users/sppatankar/Documents/MATLAB/humanStructureFunction'))

model = 'prop_prob';
data = load(fullfile(base_path, 'v3/Data', strcat(model, '_processed.mat')));


%% Load and preprocess data

n = size(data.compressibilities, 2); 

compressibilities = real(data.compressibilities);
degrees_of_freedom = data.degrees_of_freedom;
betti_dim_1 = data.betti_dim_1;
betti_dim_2 = data.betti_dim_2;
betti_dim_3 = data.betti_dim_3;

compressibilities_err = std(compressibilities, 'omitnan');
degrees_of_freedom_err = std(degrees_of_freedom, 'omitnan');
betti_dim_1_err = std(betti_dim_1, 'omitnan');
betti_dim_2_err = std(betti_dim_2, 'omitnan');
betti_dim_3_err = std(betti_dim_3, 'omitnan');

compressibilities_mean = mean(compressibilities, 'omitnan');
degrees_of_freedom_mean = mean(degrees_of_freedom, 'omitnan');
betti_dim_1_mean = mean(betti_dim_1, 'omitnan');
betti_dim_2_mean = mean(betti_dim_2, 'omitnan');
betti_dim_3_mean = mean(betti_dim_3, 'omitnan');


%% Make plots

figure;
hold on
plot(1:n, compressibilities, 'Color', [0.7, 0.7, 0.7], ...
    'LineWidth', 1);
plot(1:n, compressibilities_mean, 'LineWidth', 2, ...
    'Color', [0, 0, 0]);
hold off
xlabel('Nodes', 'FontSize', 20);
ylabel('Compressibility', 'FontSize', 20);
title('Prop. Prob.: Compressibility');
prettify

figure;
hold on
plot(1:n, degrees_of_freedom, 'Color', [0.7, 0.7, 0.7], ...
    'LineWidth', 1);
plot(1:n, degrees_of_freedom_mean, 'LineWidth', 2, ...
    'Color', [0, 0, 0]);
hold off
xlabel('Nodes', 'FontSize', 20);
ylabel('2n - e', 'FontSize', 20);
title('Prop. Prob.: Degrees of Freedom');
prettify

figure;
hold on
plot(1:n, betti_dim_1, 'Color', [0.7, 0.7, 0.7], ...
    'LineWidth', 1);
plot(1:n, betti_dim_1_mean, 'LineWidth', 2, ...
    'Color', [0, 0, 0]);
hold off
xlabel('Nodes', 'FontSize', 20);
ylabel('Dimension 1 Betti Number', 'FontSize', 20);
title('Prop. Prob.: Dimension 1 Betti Number');
prettify