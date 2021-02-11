%% Five Dimensional Curiosity Scale vs. Compressibility Peak


clc; clear; close all;

load('/Volumes/My Passport/Curiosity/v2/Data/KNOT/all_KNOT_processed.mat')
load('/Volumes/My Passport/Curiosity/Data/five_D.mat')

plot_compressibility = zeros(1, size(compressibilities, 1));
plot_betti_dim_1 = zeros(1, size(betti_dim_1, 1));
% plot_betti_dim_2 = zeros(1, size(betti_dim_2, 1));
% plot_betti_dim_3 = zeros(1, size(betti_dim_3, 1));

for i = 1:size(compressibilities, 1)
    subj_comp = compressibilities(i, :);
    idx = find(isnan(subj_comp) == 1);
    try 
        idx = idx(2) - 1; % first is always NaN, pick idx before next NaN
    catch
        idx = size(compressibilities, 2);
    end
    plot_compressibility(i) = compressibilities(i, idx);
    plot_betti_dim_1(i) = betti_dim_1(i, idx);
end

plot_compressibility = plot_compressibility';
plot_betti_dim_1 = plot_betti_dim_1';
joyous_exploration = JE_5D';
deprivation_sensitivity = DS_5D';
stress_tolerance = ST_5D';
social_curiosity = SC_5D';
thrill_seeking = TS_5D';

%% Compressibility Plots 


% 
% [r_C_JE, p_C_JE] = corr(plot_compressibility, joyous_exploration, 'type', 'Spearman');
% fprintf('r = %f, p = %f. Compressibility ~ JE\n', r_C_JE, p_C_JE);
% coeffs = polyfit(plot_compressibility, joyous_exploration, 1);
% x = linspace(min(plot_compressibility), max(plot_compressibility), 1000);
% y = polyval(coeffs, x);
% f = figure('color', 'w');
% hold on;
% scatter(plot_compressibility, joyous_exploration, 'b', 'filled');
% plot(x, y, 'LineWidth', 2);
% xlabel('Compressibility', 'FontSize', 15);
% ylabel('Joyous Exploration', 'FontSize', 15);
% title('JE vs. Compressibility', 'FontSize', 15);
% hold off;
% 
% [r_C_DS, p_C_DS] = corr(plot_compressibility, deprivation_sensitivity, 'type', 'Spearman');
% fprintf('r = %f, p = %f. Compressibility ~ DS\n', r_C_DS, p_C_DS);
% coeffs = polyfit(plot_compressibility, deprivation_sensitivity, 1);
% x = linspace(min(plot_compressibility), max(plot_compressibility), 1000);
% y = polyval(coeffs, x);
% f = figure('color', 'w');
% hold on;
% scatter(plot_compressibility, deprivation_sensitivity, 'b', 'filled');
% plot(x, y, 'LineWidth', 2);
% xlabel('Compressibility', 'FontSize', 15);
% ylabel('Deprivation Sensitivity', 'FontSize', 15);
% title('DS vs. Compressibility', 'FontSize', 15);
% hold off;
% 
% [r_C_ST, p_C_ST] = corr(plot_compressibility, stress_tolerance, 'type', 'Spearman');
% fprintf('r = %f, p = %f. Compressibility ~ ST\n', r_C_ST, p_C_ST);
% coeffs = polyfit(plot_compressibility, stress_tolerance, 1);
% x = linspace(min(plot_compressibility), max(plot_compressibility), 1000);
% y = polyval(coeffs, x);
% f = figure('color', 'w');
% hold on;
% scatter(plot_compressibility, stress_tolerance, 'b', 'filled');
% plot(x, y, 'LineWidth', 2);
% xlabel('Compressibility', 'FontSize', 15);
% ylabel('Stress Tolerance', 'FontSize', 15);
% title('ST vs. Compressibility', 'FontSize', 15);
% hold off;
% 
% [r_C_SC, p_C_SC] = corr(plot_compressibility, social_curiosity, 'type', 'Spearman');
% fprintf('r = %f, p = %f. Compressibility ~ SC\n', r_C_SC, p_C_SC);
% coeffs = polyfit(plot_compressibility, social_curiosity, 1);
% x = linspace(min(plot_compressibility), max(plot_compressibility), 1000);
% y = polyval(coeffs, x);
% f = figure('color', 'w');
% hold on;
% scatter(plot_compressibility, social_curiosity, 'b', 'filled');
% plot(x, y, 'LineWidth', 2);
% xlabel('Compressibility', 'FontSize', 15);
% ylabel('Social Curiosity', 'FontSize', 15);
% title('SC vs. Compressibility', 'FontSize', 15);
% hold off;
% 
% [r_C_TS, p_C_TS] = corr(plot_compressibility, thrill_seeking, 'type', 'Spearman');
% fprintf('r = %f, p = %f. Compressibility ~ TS\n', r_C_TS, p_C_TS);
% coeffs = polyfit(plot_compressibility, thrill_seeking, 1);
% x = linspace(min(plot_compressibility), max(plot_compressibility), 1000);
% y = polyval(coeffs, x);
% f = figure('color', 'w');
% hold on;
% scatter(plot_compressibility, thrill_seeking, 'b', 'filled');
% plot(x, y, 'LineWidth', 2);
% xlabel('Compressibility', 'FontSize', 15);
% ylabel('Thrill Seeking', 'FontSize', 15);
% title('TS vs. Compressibility', 'FontSize', 15);
% hold off;

%% Betti Number Plots 

[r_C_JE, p_C_JE] = corr(plot_betti_dim_1, joyous_exploration, 'type', 'Spearman');
fprintf('r = %f, p = %f. Betti 1 ~ JE\n', r_C_JE, p_C_JE);
coeffs = polyfit(plot_betti_dim_1, joyous_exploration, 1);
x = linspace(min(plot_betti_dim_1), max(plot_betti_dim_1), 1000);
y = polyval(coeffs, x);
f = figure('color', 'w');
hold on;
scatter(plot_betti_dim_1, joyous_exploration, 'b', 'filled');
plot(x, y, 'LineWidth', 2);
xlabel('Betti Number', 'FontSize', 15);
ylabel('Joyous Exploration', 'FontSize', 15);
title('JE vs. Betti Number', 'FontSize', 15);
hold off;

[r_C_DS, p_C_DS] = corr(plot_betti_dim_1, deprivation_sensitivity, 'type', 'Spearman');
fprintf('r = %f, p = %f. Betti 1 ~ DS\n', r_C_DS, p_C_DS);
coeffs = polyfit(plot_betti_dim_1, deprivation_sensitivity, 1);
x = linspace(min(plot_betti_dim_1), max(plot_betti_dim_1), 1000);
y = polyval(coeffs, x);
f = figure('color', 'w');
hold on;
scatter(plot_betti_dim_1, deprivation_sensitivity, 'b', 'filled');
plot(x, y, 'LineWidth', 2);
xlabel('Betti Number', 'FontSize', 15);
ylabel('Deprivation Sensitivity', 'FontSize', 15);
title('DS vs. Betti Number', 'FontSize', 15);
hold off;

[r_C_ST, p_C_ST] = corr(plot_betti_dim_1, stress_tolerance, 'type', 'Spearman');
fprintf('r = %f, p = %f. Betti 1 ~ ST\n', r_C_ST, p_C_ST);
coeffs = polyfit(plot_betti_dim_1, stress_tolerance, 1);
x = linspace(min(plot_betti_dim_1), max(plot_betti_dim_1), 1000);
y = polyval(coeffs, x);
f = figure('color', 'w');
hold on;
scatter(plot_betti_dim_1, stress_tolerance, 'b', 'filled');
plot(x, y, 'LineWidth', 2);
xlabel('Betti Number', 'FontSize', 15);
ylabel('Stress Tolerance', 'FontSize', 15);
title('ST vs. Betti Number', 'FontSize', 15);
hold off;

[r_C_SC, p_C_SC] = corr(plot_betti_dim_1, social_curiosity, 'type', 'Spearman');
fprintf('r = %f, p = %f. Betti 1 ~ SC\n', r_C_SC, p_C_SC);
coeffs = polyfit(plot_betti_dim_1, social_curiosity, 1);
x = linspace(min(plot_betti_dim_1), max(plot_betti_dim_1), 1000);
y = polyval(coeffs, x);
f = figure('color', 'w');
hold on;
scatter(plot_betti_dim_1, social_curiosity, 'b', 'filled');
plot(x, y, 'LineWidth', 2);
xlabel('Betti Number', 'FontSize', 15);
ylabel('Social Curiosity', 'FontSize', 15);
title('SC vs. Betti Number', 'FontSize', 15);
hold off;

[r_C_TS, p_C_TS] = corr(plot_betti_dim_1, thrill_seeking, 'type', 'Spearman');
fprintf('r = %f, p = %f. Betti 1 ~ TS\n', r_C_TS, p_C_TS);
coeffs = polyfit(plot_betti_dim_1, thrill_seeking, 1);
x = linspace(min(plot_betti_dim_1), max(plot_betti_dim_1), 1000);
y = polyval(coeffs, x);
f = figure('color', 'w');
hold on;
scatter(plot_betti_dim_1, thrill_seeking, 'b', 'filled');
plot(x, y, 'LineWidth', 2);
xlabel('Betti Number', 'FontSize', 15);
ylabel('Thrill Seeking', 'FontSize', 15);
title('TS vs. Betti Number', 'FontSize', 15);
hold off;