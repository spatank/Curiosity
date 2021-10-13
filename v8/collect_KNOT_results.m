clc; close all; clear;

%% Preparatory

base_path = '/Volumes/My Passport/Curiosity/';
addpath(genpath(fullfile(base_path, 'Helper')))
addpath(genpath('/Users/sppatankar/Documents/MATLAB/humanStructureFunction'))
preprocessed_data_path = ...
    fullfile(base_path, 'v8/Data/KNOT/Preprocessed/'); % for adjacency matrices
files_preprocessed = dir(fullfile(preprocessed_data_path, 'subj_*.mat'));
data_path_C = fullfile(base_path, 'v8/Data/KNOT/Processed/C/');
files_C = dir(fullfile(data_path_C, 'subj_*.mat'));
data_path_Mech = fullfile(base_path, 'v8/Data/KNOT/Processed/Mech/');
files_Mech = dir(fullfile(data_path_Mech, 'subj_*.mat'));
data_path_PH = fullfile(base_path, 'v8/Data/KNOT/Processed/PH_Matlab/');
files_PH = dir(fullfile(data_path_PH, 'subj_*.mat'));

% how large is the largest network?
max_size = 0;
for i = 1:length(files_PH)
    load(fullfile(data_path_PH, files_PH(i).name));
    if length(bettis_orig) > max_size
        max_size = length(bettis_orig);
    end
end

%% Collect and store metrics

ID = NaN(length(files_C), 1);
density = NaN(length(files_C), 1);
num_nodes = NaN(length(files_C), 1);
num_edges = NaN(length(files_C), 1);
all_C = NaN(length(files_C), max_size);
all_d = NaN(length(files_C), max_size);
all_DoF_C = NaN(length(files_C), max_size);
all_Betti_0 = NaN(length(files_C), max_size);
all_Betti_1 = NaN(length(files_C), max_size);
all_Betti_2 = NaN(length(files_C), max_size);

for i = 1:length(files_C)
    struct_1 = load(strcat(preprocessed_data_path, files_preprocessed(i).name));
    ID(i) = struct_1.subj;
    [density(i), num_nodes(i), num_edges(i)] = density_und(struct_1.G);
    struct_2 = load(strcat(data_path_C, files_C(i).name));
    assert(struct_1.subj == struct_2.subj_ID);
    all_C(i, 1:num_nodes(i)) = struct_2.C_norm;
    struct_3 = load(strcat(data_path_Mech, files_Mech(i).name));
    assert(struct_2.subj_ID == struct_3.subj_ID);
    all_d(i, 1:num_nodes(i)) = struct_3.d;
    all_DoF_C(i, 1:num_nodes(i)) = struct_3.conform;
    struct_4 = load(strcat(data_path_PH, files_PH(i).name));
    % Ideally, it would be nice to assert(struct_3.subj_ID == struct_4.subj_ID);
    all_Betti_0(i, 1:num_nodes(i)) = struct_4.bettis_orig(1, :);
    all_Betti_1(i, 1:num_nodes(i)) = struct_4.bettis_orig(2, :);
    all_Betti_2(i, 1:num_nodes(i)) = struct_4.bettis_orig(3, :);
end

ID = int64(ID); % convert to int type for downstream R `merge' usage

%% Store for analysis in R

clearvars -except ID density num_nodes num_edges ...
    all_C all_d all_DoF_C ...
    all_Betti_0 all_Betti_1 all_Betti_2

save('/Volumes/My Passport/Curiosity/v8/Data/KNOT/processed_KNOT_data.mat');
