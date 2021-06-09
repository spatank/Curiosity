clc; close all; clear;

%% Preparatory

base_path = '/Volumes/My Passport/Curiosity/';
addpath(genpath(fullfile(base_path, 'Helper')))
addpath(genpath('/Users/sppatankar/Documents/MATLAB/humanStructureFunction'))
data_path_C_DoF = fullfile(base_path, 'v6/Data/KNOT/Processed/C_DoF/');
files_C_DoF = dir(fullfile(data_path_C_DoF, 'subj_*.mat'));
data_path_Betti = fullfile(base_path, 'v6/Data/KNOT/Processed/Betti/');
files_Betti = dir(fullfile(data_path_Betti, 'subj_*.mat'));

load('/Volumes/My Passport/Curiosity/Data/five_D.mat')

%% Collect and store metrics

area_under_C = NaN(1, length(files_C_DoF));
area_under_DoF = NaN(1, length(files_C_DoF));
area_under_Betti = NaN(1, length(files_C_DoF));

max_C = NaN(1, length(files_C_DoF));
max_DoF = NaN(1, length(files_C_DoF));
max_Betti = NaN(1, length(files_C_DoF));

mean_C = NaN(1, length(files_C_DoF));
mean_DoF = NaN(1, length(files_C_DoF));
mean_Betti = NaN(1, length(files_C_DoF));

slopes_Betti = NaN(1, length(files_C_DoF));

for i = 1:length(files_C_DoF)
    load(strcat(data_path_C_DoF, files_C_DoF(i).name));
    load(strcat(data_path_Betti, files_Betti(i).name));
    area_under_C(i) = trapz(1:length(C), fillmissing(C, 'movmean', 5));
    area_under_DoF(i) = trapz(1:length(DoF), fillmissing(DoF, 'movmean', 5));
    area_under_Betti(i) = trapz(1:length(bettis_orig(2, :)), ...
        fillmissing(bettis_orig(2, :), 'movmean', 5));
    max_C(i) = max(C);
    max_DoF(i) = max(DoF);
    max_Betti(i) = max(bettis_orig(2, :));
    mean_C(i) = mean(C, 'omitnan');
    mean_DoF(i) = mean(DoF, 'omitnan');
    mean_Betti(i) = mean(bettis_orig(2, :), 'omitnan');
    coeffs = polyfit(1:length(bettis_orig(2, :)), bettis_orig(2, :), 1);
    slopes_Betti(i) = coeffs(1);
end

clearvars -except area_under_C area_under_DoF area_under_Betti max_C ...
    max_DoF max_Betti mean_C mean_DoF mean_Betti slopes_Betti ...
    JE_5D DS_5D ST_5D SC_5D TS_5D
load('/Volumes/My Passport/Curiosity/Data/five_D.mat')

%% Correlations of area under curves with 5D scale

[r_C_JE, p_C_JE] = corr(area_under_C', JE_5D', 'type', 'Spearman');
fprintf('r = %f, p = %f. AUC compressibility ~ JE\n', r_C_JE, p_C_JE);

[r_C_DS, p_C_DS] = corr(area_under_C', DS_5D', 'type', 'Spearman');
fprintf('r = %f, p = %f. AUC compressibility ~ DS\n', r_C_DS, p_C_DS);
coeffs = polyfit(area_under_C, DS_5D, 1);
x = linspace(min(area_under_C), max(area_under_C), 1000);
y = polyval(coeffs, x);
f = figure('color', 'w');
hold on;
scatter(area_under_C, DS_5D, 'b', 'filled');
plot(x, y, 'LineWidth', 2);
xlabel('Area Under Compressibility Curve', 'FontSize', 15);
ylabel('Deprivation Sensitivity', 'FontSize', 15);
hold off;

[r_C_ST, p_C_ST] = corr(area_under_C', ST_5D', 'type', 'Spearman');
fprintf('r = %f, p = %f. AUC compressibility ~ ST\n', r_C_ST, p_C_ST);

[r_C_SC, p_C_SC] = corr(area_under_C', SC_5D', 'type', 'Spearman');
fprintf('r = %f, p = %f. AUC compressibility ~ SC\n', r_C_SC, p_C_SC);

[r_C_TS, p_C_TS] = corr(area_under_C', TS_5D', 'type', 'Spearman');
fprintf('r = %f, p = %f. AUC compressibility ~ TS\n', r_C_TS, p_C_TS);

[r_Betti_JE, p_Betti_JE] = corr(area_under_Betti', JE_5D', 'type', 'Spearman');
fprintf('r = %f, p = %f. AUC Betti 1 ~ JE\n', r_Betti_JE, p_Betti_JE);

[r_Betti_DS, p_Betti_DS] = corr(area_under_Betti', DS_5D', 'type', 'Spearman');
fprintf('r = %f, p = %f. AUC Betti 1 ~ DS\n', r_Betti_DS, p_Betti_DS);
coeffs = polyfit(area_under_Betti, DS_5D, 1);
x = linspace(min(area_under_Betti), max(area_under_Betti), 1000);
y = polyval(coeffs, x);
f = figure('color', 'w');
hold on;
scatter(area_under_Betti, DS_5D, 'b', 'filled');
plot(x, y, 'LineWidth', 2);
xlabel('Area Under Betti Curve', 'FontSize', 15);
ylabel('Deprivation Sensitivity', 'FontSize', 15);
hold off;

[r_Betti_ST, p_Betti_ST] = corr(area_under_Betti', ST_5D', 'type', 'Spearman');
fprintf('r = %f, p = %f. AUC Betti 1 ~ ST\n', r_Betti_ST, p_Betti_ST);

[r_Betti_SC, p_Betti_SC] = corr(area_under_Betti', SC_5D', 'type', 'Spearman');
fprintf('r = %f, p = %f. AUC Betti 1 ~ SC\n', r_Betti_SC, p_Betti_SC);

[r_Betti_TS, p_Betti_TS] = corr(area_under_Betti', TS_5D', 'type', 'Spearman');
fprintf('r = %f, p = %f. AUC Betti 1 ~ TS\n', r_Betti_TS, p_Betti_TS);

%% Correlations of maximum curves with 5D scale

[r_C_JE, p_C_JE] = corr(max_C', JE_5D', 'type', 'Spearman');
fprintf('r = %f, p = %f. Max compressibility ~ JE\n', r_C_JE, p_C_JE);

[r_C_DS, p_C_DS] = corr(max_C', DS_5D', 'type', 'Spearman');
fprintf('r = %f, p = %f. Max compressibility ~ DS\n', r_C_DS, p_C_DS);
coeffs = polyfit(max_C, DS_5D, 1);
x = linspace(min(max_C), max(max_C), 1000);
y = polyval(coeffs, x);
f = figure('color', 'w');
hold on;
scatter(max_C, DS_5D, 'b', 'filled');
plot(x, y, 'LineWidth', 2);
xlabel('Max Compressibility', 'FontSize', 15);
ylabel('Deprivation Sensitivity', 'FontSize', 15);
hold off;

[r_C_ST, p_C_ST] = corr(max_C', ST_5D', 'type', 'Spearman');
fprintf('r = %f, p = %f. Max compressibility ~ ST\n', r_C_ST, p_C_ST);

[r_C_SC, p_C_SC] = corr(max_C', SC_5D', 'type', 'Spearman');
fprintf('r = %f, p = %f. Max compressibility ~ SC\n', r_C_SC, p_C_SC);

[r_C_TS, p_C_TS] = corr(max_C', TS_5D', 'type', 'Spearman');
fprintf('r = %f, p = %f. Max compressibility ~ TS\n', r_C_TS, p_C_TS);

[r_Betti_JE, p_Betti_JE] = corr(max_Betti', JE_5D', 'type', 'Spearman');
fprintf('r = %f, p = %f. Max Betti 1 ~ JE\n', r_Betti_JE, p_Betti_JE);

[r_Betti_DS, p_Betti_DS] = corr(max_Betti', DS_5D', 'type', 'Spearman');
fprintf('r = %f, p = %f. Max Betti 1 ~ DS\n', r_Betti_DS, p_Betti_DS);
coeffs = polyfit(max_Betti, DS_5D, 1);
x = linspace(min(max_Betti), max(max_Betti), 1000);
y = polyval(coeffs, x);
f = figure('color', 'w');
hold on;
scatter(max_Betti, DS_5D, 'b', 'filled');
plot(x, y, 'LineWidth', 2);
xlabel('Max Betti Number', 'FontSize', 15);
ylabel('Deprivation Sensitivity', 'FontSize', 15);
hold off;

[r_Betti_ST, p_Betti_ST] = corr(max_Betti', ST_5D', 'type', 'Spearman');
fprintf('r = %f, p = %f. Max Betti 1 ~ ST\n', r_Betti_ST, p_Betti_ST);

[r_Betti_SC, p_Betti_SC] = corr(max_Betti', SC_5D', 'type', 'Spearman');
fprintf('r = %f, p = %f. Max Betti 1 ~ SC\n', r_Betti_SC, p_Betti_SC);

[r_Betti_TS, p_Betti_TS] = corr(max_Betti', TS_5D', 'type', 'Spearman');
fprintf('r = %f, p = %f. Max Betti 1 ~ TS\n', r_Betti_TS, p_Betti_TS);

%% Correlations of maximum curves with 5D scale

[r_C_JE, p_C_JE] = corr(mean_C', JE_5D', 'type', 'Spearman');
fprintf('r = %f, p = %f. Mean compressibility ~ JE\n', r_C_JE, p_C_JE);

[r_C_DS, p_C_DS] = corr(mean_C', DS_5D', 'type', 'Spearman');
fprintf('r = %f, p = %f. Mean compressibility ~ DS\n', r_C_DS, p_C_DS);
coeffs = polyfit(mean_C, DS_5D, 1);
x = linspace(min(mean_C), max(mean_C), 1000);
y = polyval(coeffs, x);
f = figure('color', 'w');
hold on;
scatter(mean_C, DS_5D, 'b', 'filled');
plot(x, y, 'LineWidth', 2);
xlabel('Mean Compressibility', 'FontSize', 15);
ylabel('Deprivation Sensitivity', 'FontSize', 15);
hold off;

[r_C_ST, p_C_ST] = corr(mean_C', ST_5D', 'type', 'Spearman');
fprintf('r = %f, p = %f. Mean compressibility ~ ST\n', r_C_ST, p_C_ST);

[r_C_SC, p_C_SC] = corr(mean_C', SC_5D', 'type', 'Spearman');
fprintf('r = %f, p = %f. Mean compressibility ~ SC\n', r_C_SC, p_C_SC);

[r_C_TS, p_C_TS] = corr(mean_C', TS_5D', 'type', 'Spearman');
fprintf('r = %f, p = %f. Mean compressibility ~ TS\n', r_C_TS, p_C_TS);

[r_Betti_JE, p_Betti_JE] = corr(mean_Betti', JE_5D', 'type', 'Spearman');
fprintf('r = %f, p = %f. Mean Betti 1 ~ JE\n', r_Betti_JE, p_Betti_JE);

[r_Betti_DS, p_Betti_DS] = corr(mean_Betti', DS_5D', 'type', 'Spearman');
fprintf('r = %f, p = %f. Mean Betti 1 ~ DS\n', r_Betti_DS, p_Betti_DS);
coeffs = polyfit(mean_Betti, DS_5D, 1);
x = linspace(min(mean_Betti), max(mean_Betti), 1000);
y = polyval(coeffs, x);
f = figure('color', 'w');
hold on;
scatter(mean_Betti, DS_5D, 'b', 'filled');
plot(x, y, 'LineWidth', 2);
xlabel('Mean Betti Number', 'FontSize', 15);
ylabel('Deprivation Sensitivity', 'FontSize', 15);
hold off;

[r_Betti_ST, p_Betti_ST] = corr(mean_Betti', ST_5D', 'type', 'Spearman');
fprintf('r = %f, p = %f. Mean Betti 1 ~ ST\n', r_Betti_ST, p_Betti_ST);

[r_Betti_SC, p_Betti_SC] = corr(mean_Betti', SC_5D', 'type', 'Spearman');
fprintf('r = %f, p = %f. Mean Betti 1 ~ SC\n', r_Betti_SC, p_Betti_SC);

[r_Betti_TS, p_Betti_TS] = corr(mean_Betti', TS_5D', 'type', 'Spearman');
fprintf('r = %f, p = %f. Mean Betti 1 ~ TS\n', r_Betti_TS, p_Betti_TS);

%% Correlations of Betti curve slopes with 5D scale

[r_Betti_JE, p_Betti_JE] = corr(slopes_Betti', JE_5D', 'type', 'Spearman');
fprintf('r = %f, p = %f. Slope Betti 1 ~ JE\n', r_Betti_JE, p_Betti_JE);

[r_Betti_DS, p_Betti_DS] = corr(slopes_Betti', DS_5D', 'type', 'Spearman');
fprintf('r = %f, p = %f. Slope Betti 1 ~ DS\n', r_Betti_DS, p_Betti_DS);
coeffs = polyfit(slopes_Betti, DS_5D, 1);
x = linspace(min(slopes_Betti), max(slopes_Betti), 1000);
y = polyval(coeffs, x);
f = figure('color', 'w');
hold on;
scatter(slopes_Betti, DS_5D, 'b', 'filled');
plot(x, y, 'LineWidth', 2);
xlabel('Slope Betti Curve', 'FontSize', 15);
ylabel('Deprivation Sensitivity', 'FontSize', 15);
hold off;

[r_Betti_ST, p_Betti_ST] = corr(slopes_Betti', ST_5D', 'type', 'Spearman');
fprintf('r = %f, p = %f. Slope Betti 1 ~ ST\n', r_Betti_ST, p_Betti_ST);

[r_Betti_SC, p_Betti_SC] = corr(slopes_Betti', SC_5D', 'type', 'Spearman');
fprintf('r = %f, p = %f. Slope Betti 1 ~ SC\n', r_Betti_SC, p_Betti_SC);

[r_Betti_TS, p_Betti_TS] = corr(slopes_Betti', TS_5D', 'type', 'Spearman');
fprintf('r = %f, p = %f. Slope Betti 1 ~ TS\n', r_Betti_TS, p_Betti_TS);