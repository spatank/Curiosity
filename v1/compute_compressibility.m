clc; close all; clear;

%% Compressibility for Growing Knowledge Networks

data_path = '/Volumes/My Passport/Curiosity/Data/KNOT_data_processed';
files = dir(fullfile(data_path, 'subj_*.mat'));

setting = 7;
num_pairs = 100;

% how large is the largest network?
max_time_steps = 0;
for i = 1:length(files)
    load(fullfile(data_path, files(i).name));
    if length(all_adj) > max_time_steps
        max_time_steps = length(all_adj);
    end
end

%% Compile Compressibility Values

% initialize empty matrix of compressibility values
C_all = NaN(length(files), max_time_steps);

for i = 1:length(files)
    fprintf('Subject %d of %d.\n', i, length(files))
    load(fullfile(data_path, files(i).name));
    for j = 1:length(all_adj)
        G = full(all_adj{1, j});
        G(G > 0) = 1; % binarize
        try
            [S, S_low, clusters, Gs] = rate_distortion_upper_info(G, setting, num_pairs);
            C_all(i,j) = mean(S(end) - S);
        catch % assert that invalid return stays NaN
            C_all(i,j) = NaN;
        end
    end
end

%% Plot Time vs. Compressibility (Binary Networks)

clc; clear; close all;

load('/Volumes/My Passport/Curiosity/Data/C_all_bin.mat')

% some compressibility values are complex
C_all(imag(C_all) ~= 0) = NaN;
% some are implausibly large positive/negative numbers
C_all(abs(C_all) > 100) = NaN;

time_steps = 1:length(C_all);
C_mean = mean(C_all, 'omitnan');

f = figure('color', 'w');
scatter(time_steps, C_mean, 'b', 'filled');
xlabel('Time Steps', 'FontSize', 15);
ylabel('Compressibility', 'FontSize', 15);
title('Time vs. Compressibility (Binary Networks)', 'FontSize', 15);

%% Plot Time vs. Compressibility (Wikipedia2Vec Networks)

clc; clear; close all;

load('/Volumes/My Passport/Curiosity/Data/C_all.mat')

% some compressibility values are complex
C_all(imag(C_all) ~= 0) = NaN;
% some are implausibly large positive/negative numbers
C_all(abs(C_all) > 100) = NaN;

time_steps = 1:length(C_all);
C_mean = mean(C_all, 'omitnan');

f = figure('color', 'w');
scatter(time_steps, C_mean, 'b', 'filled');
xlabel('Time Steps', 'FontSize', 15);
ylabel('Compressibility', 'FontSize', 15);
title('Time vs. Compressibility', 'FontSize', 15);

% increase
% increase then plateau
% increase then decrease
% inverted-U
% two peaks
% three peaks

%% Five Dimensional Curiosity Scale vs. Compressibility Peak


clc; clear; close all;

load('/Volumes/My Passport/Curiosity/Data/C_all.mat')
load('/Volumes/My Passport/Curiosity/Data/five_D.mat')

C_all = real(C_all);
C_all(C_all > 5 | C_all < -5) = NaN;

plot_compressibility = nanmax(C_all, [], 2);
% plot_compressibility = mean(C_all, 2, 'omitnan');

joyous_exploration = JE_5D;
deprivation_sensitivity = DS_5D;
stress_tolerance = ST_5D;
social_curiosity = SC_5D;
thrill_seeking = TS_5D;

[r_C_JE, p_C_JE] = corr(plot_compressibility, joyous_exploration', 'type', 'Spearman');
fprintf('r = %f, p = %f. Compressibility ~ JE\n', r_C_JE, p_C_JE);
coeffs = polyfit(plot_compressibility', joyous_exploration, 1);
x = linspace(min(plot_compressibility), max(plot_compressibility), 1000);
y = polyval(coeffs, x);
f = figure('color', 'w');
hold on;
scatter(plot_compressibility, joyous_exploration, 'b', 'filled');
plot(x, y, 'LineWidth', 2);
xlabel('Compressibility', 'FontSize', 15);
ylabel('Joyous Exploration', 'FontSize', 15);
title('JE vs. Compressibility', 'FontSize', 15);
hold off;

[r_C_DS, p_C_DS] = corr(plot_compressibility, deprivation_sensitivity', 'type', 'Spearman');
fprintf('r = %f, p = %f. Compressibility ~ DS\n', r_C_DS, p_C_DS);
coeffs = polyfit(plot_compressibility', deprivation_sensitivity, 1);
x = linspace(min(plot_compressibility), max(plot_compressibility), 1000);
y = polyval(coeffs, x);
f = figure('color', 'w');
hold on;
scatter(plot_compressibility, deprivation_sensitivity, 'b', 'filled');
plot(x, y, 'LineWidth', 2);
xlabel('Compressibility', 'FontSize', 15);
ylabel('Deprivation Sensitivity', 'FontSize', 15);
title('DS vs. Compressibility', 'FontSize', 15);
hold off;

[r_C_ST, p_C_ST] = corr(plot_compressibility, stress_tolerance', 'type', 'Spearman');
fprintf('r = %f, p = %f. Compressibility ~ ST\n', r_C_ST, p_C_ST);
coeffs = polyfit(plot_compressibility', stress_tolerance, 1);
x = linspace(min(plot_compressibility), max(plot_compressibility), 1000);
y = polyval(coeffs, x);
f = figure('color', 'w');
hold on;
scatter(plot_compressibility, stress_tolerance, 'b', 'filled');
plot(x, y, 'LineWidth', 2);
xlabel('Compressibility', 'FontSize', 15);
ylabel('Stress Tolerance', 'FontSize', 15);
title('ST vs. Compressibility', 'FontSize', 15);
hold off;

[r_C_SC, p_C_SC] = corr(plot_compressibility, social_curiosity', 'type', 'Spearman');
fprintf('r = %f, p = %f. Compressibility ~ SC\n', r_C_SC, p_C_SC);
coeffs = polyfit(plot_compressibility', social_curiosity, 1);
x = linspace(min(plot_compressibility), max(plot_compressibility), 1000);
y = polyval(coeffs, x);
f = figure('color', 'w');
hold on;
scatter(plot_compressibility, social_curiosity, 'b', 'filled');
plot(x, y, 'LineWidth', 2);
xlabel('Compressibility', 'FontSize', 15);
ylabel('Social Curiosity', 'FontSize', 15);
title('SC vs. Compressibility', 'FontSize', 15);
hold off;

[r_C_TS, p_C_TS] = corr(plot_compressibility, thrill_seeking', 'type', 'Spearman');
fprintf('r = %f, p = %f. Compressibility ~ TS\n', r_C_TS, p_C_TS);
coeffs = polyfit(plot_compressibility', thrill_seeking, 1);
x = linspace(min(plot_compressibility), max(plot_compressibility), 1000);
y = polyval(coeffs, x);
f = figure('color', 'w');
hold on;
scatter(plot_compressibility, thrill_seeking, 'b', 'filled');
plot(x, y, 'LineWidth', 2);
xlabel('Compressibility', 'FontSize', 15);
ylabel('Thrill Seeking', 'FontSize', 15);
title('TS vs. Compressibility', 'FontSize', 15);
hold off;
