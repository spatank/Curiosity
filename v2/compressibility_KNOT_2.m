clc; close all; clear;

%% Load all KNOT time series

load('/Volumes/My Passport/Curiosity/v2/Data/KNOT/all_KNOT_processed.mat')
num_subjects = size(compressibilities, 1);
all_lengths = zeros(1, num_subjects);

%% Re-sample compressibility curves

compressibilities = compressibilities(:, 2:end);
betti_dim_1 = betti_dim_1(:, 2:end);
betti_dim_2 = betti_dim_2(:, 2:end);
betti_dim_3 = betti_dim_3(:, 2:end);

for i = 1:num_subjects
    subj_betti_dim_1 = betti_dim_1(i, :);
    idx = find(isnan(subj_betti_dim_1) == 1);
    try
        all_lengths(i) = idx(1) - 1;
    catch
        fprintf('Largest network for subject %d.\n', i);
        all_lengths(i) = length(subj_betti_dim_1);
    end
end

new_length = round(mean(all_lengths));
resampled_comp = zeros(num_subjects, new_length);
resampled_betti_dim_1 = zeros(num_subjects, new_length);
resampled_betti_dim_2 = zeros(num_subjects, new_length);
resampled_betti_dim_3 = zeros(num_subjects, new_length);

for i = 1:num_subjects
    subj_comp = compressibilities(i, 1:all_lengths(i));
    subj_betti_dim_1 = betti_dim_1(i, 1:all_lengths(i));
    subj_betti_dim_2 = betti_dim_2(i, 1:all_lengths(i));
    subj_betti_dim_3 = betti_dim_3(i, 1:all_lengths(i));
    if any(isnan(subj_comp))
        fprintf('NaN value for subject %d.\n', i);
    end
    subj_len = all_lengths(i);
    resampled_comp(i, :) = interp1(1:length(subj_comp), subj_comp, ...
        linspace(1, length(subj_comp), new_length), 'pchip');
    resampled_betti_dim_1(i, :) = interp1(1:length(subj_betti_dim_1), subj_betti_dim_1, ...
        linspace(1, length(subj_comp), new_length), 'pchip');
    resampled_betti_dim_2(i, :) = interp1(1:length(subj_betti_dim_2), subj_betti_dim_2, ...
        linspace(1, length(subj_comp), new_length), 'pchip');
    resampled_betti_dim_3(i, :) = interp1(1:length(subj_betti_dim_3), subj_betti_dim_3, ...
        linspace(1, length(subj_comp), new_length), 'pchip');
end

%% Plot Time vs. Compressibility 

new_length = size(resampled_comp, 2); 

nodes = 1:new_length;
compressibilities_raw = mean(resampled_comp, 'omitnan');

betti_dim_1_raw = mean(resampled_betti_dim_1, 'omitnan');

betti_dim_2_raw = mean(resampled_betti_dim_2, 'omitnan');

betti_dim_3_raw = mean(resampled_betti_dim_3, 'omitnan');



% FIGURE OF INTEREST
figure;
plot(compressibilities_raw, betti_dim_1_raw, 'LineWidth', 2, ...
    'Color', [0, 0, 0]);
ylim([-1, 6])
xlabel('Compressibility', 'FontSize', 20);
ylabel('Betti Number (Dimension 1)', 'FontSize', 20);
title('KNOT', 'FontSize', 15);
prettify



figure;
plot(1:new_length, compressibilities_raw, 'LineWidth', 2, ...
    'Color', [0, 0, 0]);
xlabel('Nodes', 'FontSize', 20);
ylabel('Compressibility', 'FontSize', 20);
title('KNOT');
prettify

figure;
plot(1:new_length, betti_dim_1_raw, 'LineWidth', 2, ...
    'Color', [0, 0, 0]);

xlabel('Nodes', 'FontSize', 20);
ylabel('Betti Number (Dimension 1)', 'FontSize', 20);
title('KNOT');
prettify
