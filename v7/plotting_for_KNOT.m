clc; clear;

%% Individual

base_path = '/Volumes/My Passport/Curiosity/';
addpath(genpath(fullfile(base_path, 'Helper')))
data_path = '/Volumes/My Passport/Curiosity/v7/Data/KNOT/';
subj_ID = 363;
load(strcat(data_path, 'Processed/C_clust_coef/subj_', string(subj_ID), '_C_clust_coef.mat'))
load(strcat(data_path, 'Processed/Betti/subj_', string(subj_ID), '_bettis.mat'))
load(strcat(data_path, 'Processed/Mech/subj_', string(subj_ID), '_mech.mat'))
n = size(conform, 2);
num_iters = size(conform_edge_rewired, 1);

figure;
hold on
plot(1:n, C, 'Color', [0, 0, 0], ...
    'LineWidth', 1);
errorbar(1:n, mean(C_edge_rewired, 'omitnan'), std(C_edge_rewired, 'omitnan')/num_iters, ...
    'MarkerSize', 2, 'LineWidth', 1, 'Color', [0.8500, 0.3250, 0.0980], ...
    'MarkerEdgeColor',[0.8500, 0.3250, 0.0980], ...
    'MarkerFaceColor',[0.8500, 0.3250, 0.0980]);
errorbar(1:n, mean(C_latticized, 'omitnan'), std(C_latticized, 'omitnan')/num_iters, ...
    'MarkerSize', 2, 'LineWidth', 1, 'Color', [0.9290, 0.6940, 0.1250], ...
    'MarkerEdgeColor',[0.9290, 0.6940, 0.1250], ...
    'MarkerFaceColor',[0.9290, 0.6940, 0.1250]);
hold off
xlabel('Nodes', 'FontSize', 20);
ylabel('Compressibility', 'FontSize', 20);
legend('Original', 'Edges Rewired', 'Latticized', ...
    'Location', 'NW');
title(strcat({'Subj. '}, string(subj_ID), {': Compressibility'}));
prettify

figure;
hold on
plot(1:n, clust_coef, 'Color', [0, 0, 0], ...
    'LineWidth', 1);
errorbar(1:n, mean(clust_coef_edge_rewired, 'omitnan'), std(clust_coef_edge_rewired, 'omitnan')/num_iters, ...
    'MarkerSize', 2, 'LineWidth', 1, 'Color', [0.8500, 0.3250, 0.0980], ...
    'MarkerEdgeColor',[0.8500, 0.3250, 0.0980], ...
    'MarkerFaceColor',[0.8500, 0.3250, 0.0980]);
errorbar(1:n, mean(clust_coef_latticized, 'omitnan'), std(clust_coef_latticized, 'omitnan')/num_iters, ...
    'MarkerSize', 2, 'LineWidth', 1, 'Color', [0.9290, 0.6940, 0.1250], ...
    'MarkerEdgeColor',[0.9290, 0.6940, 0.1250], ...
    'MarkerFaceColor',[0.9290, 0.6940, 0.1250]);
hold off
xlabel('Nodes', 'FontSize', 20);
ylabel('Mean Clustering Coefficient', 'FontSize', 20);
legend('Original', 'Edges Rewired', 'Latticized', ...
    'Location', 'NW');
title(strcat({'Subj. '}, string(subj_ID), {': Clustering Coefficient'}));
prettify

figure;
hold on
plot(1:n, bettis_orig(2, :), 'Color', [0, 0, 0], ...
    'LineWidth', 1);
errorbar(1:n, mean(bettis_1_edges_rewired, 'omitnan'), std(bettis_1_edges_rewired, 'omitnan')/num_iters, ...
    'MarkerSize', 2, 'LineWidth', 1, 'Color', [0.8500, 0.3250, 0.0980], ...
    'MarkerEdgeColor',[0.8500, 0.3250, 0.0980], ...
    'MarkerFaceColor',[0.8500, 0.3250, 0.0980]);
errorbar(1:n, mean(bettis_1_latticized, 'omitnan'), std(bettis_1_latticized, 'omitnan')/num_iters, ...
    'MarkerSize', 2, 'LineWidth', 1, 'Color', [0.9290, 0.6940, 0.1250], ...
    'MarkerEdgeColor',[0.9290, 0.6940, 0.1250], ...
    'MarkerFaceColor',[0.9290, 0.6940, 0.1250]);
hold off
xlabel('Nodes', 'FontSize', 20);
ylabel('Betti Number', 'FontSize', 20);
legend('Original', 'Edges Rewired', 'Latticized', ...
    'Location', 'NW');
title(strcat({'Subj. '}, string(subj_ID), {': Betti Curve'}));
prettify

figure;
hold on
plot(1:n, d, 'Color', [0, 0, 0], ...
    'LineWidth', 1);
errorbar(1:n, mean(d_edge_rewired, 'omitnan'), std(DoF_edge_rewired, 'omitnan')/num_iters, ...
    'MarkerSize', 2, 'LineWidth', 1, 'Color', [0.8500, 0.3250, 0.0980], ...
    'MarkerEdgeColor',[0.8500, 0.3250, 0.0980], ...
    'MarkerFaceColor',[0.8500, 0.3250, 0.0980]);
errorbar(1:n, mean(d_latticized, 'omitnan'), std(DoF_latticized, 'omitnan')/num_iters, ...
    'MarkerSize', 2, 'LineWidth', 1, 'Color', [0.9290, 0.6940, 0.1250], ...
    'MarkerEdgeColor',[0.9290, 0.6940, 0.1250], ...
    'MarkerFaceColor',[0.9290, 0.6940, 0.1250]);
hold off
xlabel('Nodes', 'FontSize', 20);
ylabel('Embedding Dimensionality', 'FontSize', 20);
legend('Original', 'Edges Rewired', 'Latticized', ...
    'Location', 'NW');
title(strcat({'Subj. '}, string(subj_ID), {': d'}));
prettify

figure;
hold on
plot(1:n, DoF, 'Color', [0, 0, 0], ...
    'LineWidth', 1);
errorbar(1:n, mean(DoF_edge_rewired, 'omitnan'), std(DoF_edge_rewired, 'omitnan')/num_iters, ...
    'MarkerSize', 2, 'LineWidth', 1, 'Color', [0.8500, 0.3250, 0.0980], ...
    'MarkerEdgeColor',[0.8500, 0.3250, 0.0980], ...
    'MarkerFaceColor',[0.8500, 0.3250, 0.0980]);
errorbar(1:n, mean(DoF_latticized, 'omitnan'), std(DoF_latticized, 'omitnan')/num_iters, ...
    'MarkerSize', 2, 'LineWidth', 1, 'Color', [0.9290, 0.6940, 0.1250], ...
    'MarkerEdgeColor',[0.9290, 0.6940, 0.1250], ...
    'MarkerFaceColor',[0.9290, 0.6940, 0.1250]);
hold off
xlabel('Nodes', 'FontSize', 20);
ylabel('DoF', 'FontSize', 20);
legend('Original', 'Edges Rewired', 'Latticized', ...
    'Location', 'NW');
title(strcat({'Subj. '}, string(subj_ID), {': DoF'}));
prettify

figure;
hold on
plot(1:n, rigid, 'Color', [0, 0, 0], ...
    'LineWidth', 1);
errorbar(1:n, mean(rigid_edge_rewired, 'omitnan'), std(DoF_edge_rewired, 'omitnan')/num_iters, ...
    'MarkerSize', 2, 'LineWidth', 1, 'Color', [0.8500, 0.3250, 0.0980], ...
    'MarkerEdgeColor',[0.8500, 0.3250, 0.0980], ...
    'MarkerFaceColor',[0.8500, 0.3250, 0.0980]);
errorbar(1:n, mean(rigid_latticized, 'omitnan'), std(DoF_latticized, 'omitnan')/num_iters, ...
    'MarkerSize', 2, 'LineWidth', 1, 'Color', [0.9290, 0.6940, 0.1250], ...
    'MarkerEdgeColor',[0.9290, 0.6940, 0.1250], ...
    'MarkerFaceColor',[0.9290, 0.6940, 0.1250]);
hold off
xlabel('Nodes', 'FontSize', 20);
ylabel('Rigid Body DoF', 'FontSize', 20);
legend('Original', 'Edges Rewired', 'Latticized', ...
    'Location', 'NW');
title(strcat({'Subj. '}, string(subj_ID), {': Rigid Body DoF'}));
prettify

figure;
hold on
plot(1:n, conform, 'Color', [0, 0, 0], ...
    'LineWidth', 1);
errorbar(1:n, mean(conform_edge_rewired, 'omitnan'), std(DoF_edge_rewired, 'omitnan')/num_iters, ...
    'MarkerSize', 2, 'LineWidth', 1, 'Color', [0.8500, 0.3250, 0.0980], ...
    'MarkerEdgeColor',[0.8500, 0.3250, 0.0980], ...
    'MarkerFaceColor',[0.8500, 0.3250, 0.0980]);
errorbar(1:n, mean(conform_latticized, 'omitnan'), std(DoF_latticized, 'omitnan')/num_iters, ...
    'MarkerSize', 2, 'LineWidth', 1, 'Color', [0.9290, 0.6940, 0.1250], ...
    'MarkerEdgeColor',[0.9290, 0.6940, 0.1250], ...
    'MarkerFaceColor',[0.9290, 0.6940, 0.1250]);
hold off
xlabel('Nodes', 'FontSize', 20);
ylabel('Conformable DoF', 'FontSize', 20);
legend('Original', 'Edges Rewired', 'Latticized', ...
    'Location', 'NW');
title(strcat({'Subj. '}, string(subj_ID), {': Conformable DoF'}));
prettify

%% All

clc; clear;

base_path = '/Volumes/My Passport/Curiosity/';
addpath(genpath(fullfile(base_path, 'Helper')))
addpath(genpath('/Users/sppatankar/Documents/MATLAB/humanStructureFunction'))
data_path = fullfile(base_path, 'v7/Data/KNOT/Preprocessed/');
files = dir(fullfile(data_path, 'subj_*.mat'));

num_subjs = length(files);

% how large is the largest network?
max_size = 0;
for i = 1:length(files)
    load(fullfile(data_path, files(i).name));
    if size(G, 1) > max_size
        max_size = size(G, 1);
    end
end

proc_data_path = '/Volumes/My Passport/Curiosity/v7/Data/KNOT/Processed';
num_iters = 25;

C_clust_coef_files = dir(fullfile(proc_data_path, 'C_clust_coef', 'subj_*.mat'));

all_clust_coef = NaN(length(files), max_size);
all_clust_coef_rewired = NaN(length(files), max_size);
all_clust_coef_latticized = NaN(length(files), max_size);

for i = 1:length(C_clust_coef_files)
    load(fullfile(proc_data_path, 'C_clust_coef', C_clust_coef_files(i).name));
    all_clust_coef(i, 1:length(clust_coef)) = clust_coef;
    all_clust_coef_rewired(i, 1:length(clust_coef)) = mean(clust_coef_edge_rewired);
    all_clust_coef_latticized(i, 1:length(clust_coef)) = mean(clust_coef_latticized);
end

figure;
hold on
plot(1:max_size, all_clust_coef, 'Color', [0.7, 0.7, 0.7], ...
    'LineWidth', 1);
plot(1:max_size, mean(all_clust_coef, 'omitnan'), 'LineWidth', 2, ...
    'Color', [0, 0, 0]);
hold off
xlabel('Nodes', 'FontSize', 20);
ylabel('Mean Clustering Coefficient', 'FontSize', 20);
prettify

figure;
hold on
% plot(1:max_size, mean(all_clust_coef, 'omitnan'), 'Color', [0, 0, 0], ...
%     'LineWidth', 2);
errorbar(1:max_size, mean(all_clust_coef, 'omitnan'), std(all_clust_coef, 'omitnan')/num_subjs, ...
    'MarkerSize', 2, 'LineWidth', 1, 'Color', [0, 0, 0], ...
    'MarkerEdgeColor', [0, 0, 0], ...
    'MarkerFaceColor', [0, 0, 0]);
errorbar(1:max_size, mean(all_clust_coef_rewired, 'omitnan'), std(all_clust_coef_rewired, 'omitnan')/num_subjs, ...
    'MarkerSize', 2, 'LineWidth', 1, 'Color', [0.8500, 0.3250, 0.0980], ...
    'MarkerEdgeColor',[0.8500, 0.3250, 0.0980], ...
    'MarkerFaceColor',[0.8500, 0.3250, 0.0980]);
errorbar(1:max_size, mean(all_clust_coef_latticized, 'omitnan'), std(all_clust_coef_latticized, 'omitnan')/num_subjs, ...
    'MarkerSize', 2, 'LineWidth', 1, 'Color', [0.9290, 0.6940, 0.1250], ...
    'MarkerEdgeColor',[0.9290, 0.6940, 0.1250], ...
    'MarkerFaceColor',[0.9290, 0.6940, 0.1250]);
hold off
xlabel('Nodes', 'FontSize', 20);
ylabel('Mean Clustering Coefficient', 'FontSize', 20);
legend('Original', 'Edges Rewired', 'Latticized', ...
    'Location', 'NW');
prettify

% all_C = NaN(length(files), max_size);
% all_C_rewired = NaN(length(files), max_size);
% all_C_latticized = NaN(length(files), max_size);
% 
% for i = 1:length(C_clust_coef_files)
%     load(fullfile(proc_data_path, 'C_clust_coef', C_clust_coef_files(i).name));
%     all_C(i, 1:length(C)) = C;
%     all_C_rewired(i, 1:length(C)) = mean(C_edge_rewired);
%     all_C_latticized(i, 1:length(C)) = mean(C_latticized);
% end
% 
% figure;
% hold on
% plot(1:max_size, all_C, 'Color', [0.7, 0.7, 0.7], ...
%     'LineWidth', 1);
% plot(1:max_size, mean(all_C, 'omitnan'), 'LineWidth', 2, ...
%     'Color', [0, 0, 0]);
% hold off
% xlabel('Nodes', 'FontSize', 20);
% ylabel('Compressibility', 'FontSize', 20);
% prettify
% 
% figure;
% hold on
% % plot(1:max_size, mean(all_C, 'omitnan'), 'Color', [0, 0, 0], ...
% %     'LineWidth', 2);
% errorbar(1:max_size, mean(all_C, 'omitnan'), std(all_C, 'omitnan')/num_subjs, ...
%     'MarkerSize', 2, 'LineWidth', 1, 'Color', [0, 0, 0], ...
%     'MarkerEdgeColor', [0, 0, 0], ...
%     'MarkerFaceColor', [0, 0, 0]);
% errorbar(1:max_size, mean(all_C_rewired, 'omitnan'), std(all_C_rewired, 'omitnan')/num_subjs, ...
%     'MarkerSize', 2, 'LineWidth', 1, 'Color', [0.8500, 0.3250, 0.0980], ...
%     'MarkerEdgeColor',[0.8500, 0.3250, 0.0980], ...
%     'MarkerFaceColor',[0.8500, 0.3250, 0.0980]);
% errorbar(1:max_size, mean(all_C_latticized, 'omitnan'), std(all_C_latticized, 'omitnan')/num_subjs, ...
%     'MarkerSize', 2, 'LineWidth', 1, 'Color', [0.9290, 0.6940, 0.1250], ...
%     'MarkerEdgeColor',[0.9290, 0.6940, 0.1250], ...
%     'MarkerFaceColor',[0.9290, 0.6940, 0.1250]);
% hold off
% xlabel('Nodes', 'FontSize', 20);
% ylabel('Compressibility', 'FontSize', 20);
% legend('Original', 'Edges Rewired', 'Latticized', ...
%     'Location', 'NW');
% prettify

Betti_files = dir(fullfile(proc_data_path, 'Betti', 'subj_*.mat'));

all_Betti = NaN(length(Betti_files), max_size);
all_Betti_rewired = NaN(length(Betti_files), max_size);
all_Betti_latticized = NaN(length(Betti_files), max_size);

for i = 1:length(Betti_files)
    load(fullfile(proc_data_path, 'Betti', Betti_files(i).name));
    all_Betti(i, 1:length(bettis_orig(2, :))) = bettis_orig(2, :);
    all_Betti_rewired(i, 1:length(bettis_orig(2, :))) = ...
        mean(bettis_1_edges_rewired);
    all_Betti_latticized(i, 1:length(bettis_orig(2, :))) = ...
        mean(bettis_1_latticized);
end

figure;
hold on
plot(1:max_size, all_Betti, 'Color', [0.7, 0.7, 0.7], ...
    'LineWidth', 1);
plot(1:max_size, mean(all_Betti, 'omitnan'), 'LineWidth', 2, ...
    'Color', [0, 0, 0]);
hold off
xlabel('Nodes', 'FontSize', 20);
ylabel('Cycles', 'FontSize', 20);
prettify

figure;
hold on
% plot(1:max_size, mean(all_Betti, 'omitnan'), 'Color', [0, 0, 0], ...
%     'LineWidth', 2);
errorbar(1:max_size, mean(all_Betti, 'omitnan'), std(all_Betti, 'omitnan')/num_subjs, ...
    'MarkerSize', 2, 'LineWidth', 1, 'Color', [0, 0, 0], ...
    'MarkerEdgeColor', [0, 0, 0], ...
    'MarkerFaceColor', [0, 0, 0]);
errorbar(1:max_size, mean(all_Betti_rewired, 'omitnan'), std(all_Betti_rewired, 'omitnan')/num_subjs, ...
    'MarkerSize', 2, 'LineWidth', 1, 'Color', [0.8500, 0.3250, 0.0980], ...
    'MarkerEdgeColor', [0.8500, 0.3250, 0.0980], ...
    'MarkerFaceColor', [0.8500, 0.3250, 0.0980]);
errorbar(1:max_size, mean(all_Betti_latticized, 'omitnan'), std(all_Betti_latticized, 'omitnan')/num_subjs, ...
    'MarkerSize', 2, 'LineWidth', 1, 'Color', [0.9290, 0.6940, 0.1250], ...
    'MarkerEdgeColor', [0.9290, 0.6940, 0.1250], ...
    'MarkerFaceColor', [0.9290, 0.6940, 0.1250]);
hold off
xlabel('Nodes', 'FontSize', 20);
ylabel('Cycles', 'FontSize', 20);
legend('Original', 'Edges Rewired', 'Latticized', ...
    'Location', 'NW');
prettify

mech_files = dir(fullfile(proc_data_path, 'Mech', 'subj_*.mat'));

all_d = NaN(length(mech_files), max_size);
all_d_rewired = NaN(length(mech_files), max_size);
all_d_latticized = NaN(length(mech_files), max_size);

for i = 1:length(mech_files)
    load(fullfile(proc_data_path, 'Mech', mech_files(i).name));
    all_d(i, 1:length(d)) = d;
    all_d_rewired(i, 1:length(d)) = ...
        mean(d_edge_rewired);
    all_d_latticized(i, 1:length(d)) = ...
        mean(d_latticized);
end

all_DoF = NaN(length(mech_files), max_size);
all_DoF_rewired = NaN(length(mech_files), max_size);
all_DoF_latticized = NaN(length(mech_files), max_size);

for i = 1:length(mech_files)
    load(fullfile(proc_data_path, 'Mech', mech_files(i).name));
    all_DoF(i, 1:length(DoF)) = DoF;
    all_DoF_rewired(i, 1:length(DoF)) = ...
        mean(DoF_edge_rewired);
    all_DoF_latticized(i, 1:length(DoF)) = ...
        mean(DoF_latticized);
end

all_rigid = NaN(length(mech_files), max_size);
all_rigid_rewired = NaN(length(mech_files), max_size);
all_rigid_latticized = NaN(length(mech_files), max_size);

for i = 1:length(mech_files)
    load(fullfile(proc_data_path, 'Mech', mech_files(i).name));
    all_rigid(i, 1:length(rigid)) = rigid;
    all_rigid_rewired(i, 1:length(rigid)) = ...
        mean(rigid_edge_rewired);
    all_rigid_latticized(i, 1:length(rigid)) = ...
        mean(rigid_latticized);
end

all_conform = NaN(length(mech_files), max_size);
all_conform_rewired = NaN(length(mech_files), max_size);
all_conform_latticized = NaN(length(mech_files), max_size);

for i = 1:length(mech_files)
    load(fullfile(proc_data_path, 'Mech', mech_files(i).name));
    all_conform(i, 1:length(conform)) = conform;
    all_conform_rewired(i, 1:length(conform)) = ...
        mean(conform_edge_rewired);
    all_conform_latticized(i, 1:length(conform)) = ...
        mean(conform_latticized);
end

figure;
hold on
plot(1:max_size, all_d, 'Color', [0.7, 0.7, 0.7], ...
    'LineWidth', 1);
plot(1:max_size, (mean(all_d, 'omitnan')), 'LineWidth', 2, ...
    'Color', [0, 0, 0]);
hold off
xlabel('Nodes', 'FontSize', 20);
ylabel('Embedding Dimensionality', 'FontSize', 20);
prettify

figure;
hold on
plot(1:max_size, (mean(all_d, 'omitnan')), 'LineWidth', 2, ...
    'Color', [0, 0, 0]);
plot(1:max_size, (mean(all_d_rewired, 'omitnan')), 'LineWidth', 2, ...
    'Color', [0.8500, 0.3250, 0.0980]);
plot(1:max_size, (mean(all_d_latticized, 'omitnan')), 'LineWidth', 2, ...
    'Color', [0.9290, 0.6940, 0.1250]);
hold off
xlabel('Nodes', 'FontSize', 20);
ylabel('Embedding Dimensionality', 'FontSize', 20);
legend('Original', 'Edges Rewired', 'Latticized', ...
    'Location', 'NW');
prettify

figure;
hold on
plot(1:max_size, all_DoF, 'Color', [0.7, 0.7, 0.7], ...
    'LineWidth', 1);
plot(1:max_size, (mean(all_DoF, 'omitnan')), 'LineWidth', 2, ...
    'Color', [0, 0, 0]);
hold off
xlabel('Nodes', 'FontSize', 20);
ylabel('DoF', 'FontSize', 20);
prettify

figure;
hold on
plot(1:max_size, (mean(all_DoF, 'omitnan')), 'LineWidth', 2, ...
    'Color', [0, 0, 0]);
plot(1:max_size, (mean(all_DoF_rewired, 'omitnan')), 'LineWidth', 2, ...
    'Color', [0.8500, 0.3250, 0.0980]);
plot(1:max_size, (mean(all_DoF_latticized, 'omitnan')), 'LineWidth', 2, ...
    'Color', [0.9290, 0.6940, 0.1250]);
hold off
xlabel('Nodes', 'FontSize', 20);
ylabel('DoF', 'FontSize', 20);
legend('Original', 'Edges Rewired', 'Latticized', ...
    'Location', 'NW');
prettify

figure;
hold on
plot(1:max_size, all_rigid, 'Color', [0.7, 0.7, 0.7], ...
    'LineWidth', 1);
plot(1:max_size, (mean(all_rigid, 'omitnan')), 'LineWidth', 2, ...
    'Color', [0, 0, 0]);
hold off
xlabel('Nodes', 'FontSize', 20);
ylabel('Rigid Body DoF', 'FontSize', 20);
prettify

figure;
hold on
plot(1:max_size, (mean(all_rigid, 'omitnan')), 'LineWidth', 2, ...
    'Color', [0, 0, 0]);
plot(1:max_size, (mean(all_rigid_rewired, 'omitnan')), 'LineWidth', 2, ...
    'Color', [0.8500, 0.3250, 0.0980]);
plot(1:max_size, (mean(all_rigid_latticized, 'omitnan')), 'LineWidth', 2, ...
    'Color', [0.9290, 0.6940, 0.1250]);
hold off
xlabel('Nodes', 'FontSize', 20);
ylabel('Rigid Body DoF', 'FontSize', 20);
legend('Original', 'Edges Rewired', 'Latticized', ...
    'Location', 'NW');
prettify

figure;
hold on
plot(1:max_size, all_conform, 'Color', [0.7, 0.7, 0.7], ...
    'LineWidth', 1);
plot(1:max_size, (mean(all_conform, 'omitnan')), 'LineWidth', 2, ...
    'Color', [0, 0, 0]);
hold off
xlabel('Nodes', 'FontSize', 20);
ylabel('Conformable DoF', 'FontSize', 20);
prettify

figure;
hold on
plot(1:max_size, (mean(all_conform, 'omitnan')), 'LineWidth', 2, ...
    'Color', [0, 0, 0]);
plot(1:max_size, (mean(all_conform_rewired, 'omitnan')), 'LineWidth', 2, ...
    'Color', [0.8500, 0.3250, 0.0980]);
plot(1:max_size, (mean(all_conform_latticized, 'omitnan')), 'LineWidth', 2, ...
    'Color', [0.9290, 0.6940, 0.1250]);
hold off
xlabel('Nodes', 'FontSize', 20);
ylabel('Conformable DoF', 'FontSize', 20);
legend('Original', 'Edges Rewired', 'Latticized', ...
    'Location', 'NW');
prettify