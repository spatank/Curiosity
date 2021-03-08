clc; close all; clear;

%% Load all KNOT time series

base_path = '/Volumes/My Passport/Curiosity/';
addpath(genpath(fullfile(base_path, 'Helper')))

load('/Volumes/My Passport/Curiosity/v3/test_1.mat')
num_subjects = size(compressibilities, 1);
all_lengths = zeros(1, num_subjects);

zscore_nan = @(x) bsxfun(@rdivide, bsxfun(@minus, x, mean(x, 'omitnan')), ...
    std(x, 'omitnan'));

figure;
hold on
plot(1:length(compressibilities_raw), zscore_nan(compressibilities_raw), 'LineWidth', 2, ...
    'Color', [0, 0, 0]);
plot(1:length(degrees_of_freedom_raw), zscore(degrees_of_freedom_raw), 'LineWidth', 2, ...
    'Color', [0.8500, 0.3250, 0.0980]);
xlabel('Node', 'FontSize', 20);
ylabel('Z-Score', 'FontSize', 20);
title('KNOT', 'FontSize', 20);
legend('Compressibility', '2n - e', 'Location', 'NorthWest');
prettify
