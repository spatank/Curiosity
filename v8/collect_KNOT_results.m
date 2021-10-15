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
num_change_points = NaN(length(files_C), 1);
all_DoF_C = NaN(length(files_C), max_size);
all_Betti_0 = NaN(length(files_C), max_size);
closed_dim_0 = NaN(length(files_C), 1);
all_Betti_1 = NaN(length(files_C), max_size);
closed_dim_1 = NaN(length(files_C), 1);
all_Betti_2 = NaN(length(files_C), max_size);
closed_dim_2 = NaN(length(files_C), 1);

for i = 1:length(files_C)
    struct_1 = load(strcat(preprocessed_data_path, files_preprocessed(i).name));
    [density(i), num_nodes(i), num_edges(i)] = density_und(struct_1.G);
    struct_2 = load(strcat(data_path_C, files_C(i).name));
    assert(struct_1.subj == struct_2.subj_ID);
    all_C(i, 1:num_nodes(i)) = struct_2.C_norm;
    struct_3 = load(strcat(data_path_Mech, files_Mech(i).name));
    assert(struct_2.subj_ID == struct_3.subj_ID);
    all_d(i, 1:num_nodes(i)) = struct_3.d;
    num_change_points(i) = length(find(diff(struct_3.d) == 1));
    all_DoF_C(i, 1:num_nodes(i)) = struct_3.conform;
    struct_4 = load(strcat(data_path_PH, files_PH(i).name));
    % Ideally, it would be nice to easily assert(struct_3.subj_ID == struct_4.subj_ID);
    subj_ID = strsplit(files_PH(i).name, '_');
    subj_ID = str2double(cell2mat(subj_ID(2))); % 1: `subj', 2: ID, 3: 'PH.mat'
    assert(struct_3.subj_ID == subj_ID)
    all_Betti_0(i, 1:num_nodes(i)) = struct_4.bettis_orig(1, :);
    bars_dim_0 = struct_4.bars_orig{1, 2};
    if ~isempty(bars_dim_0)
        closed_dim_0(i) = numel(find(bars_dim_0(:, 2) ~= Inf));
    end
    all_Betti_1(i, 1:num_nodes(i)) = struct_4.bettis_orig(2, :);
    bars_dim_1 = struct_4.bars_orig{2, 2};
    if ~isempty(bars_dim_1)
        closed_dim_1(i) = numel(find(bars_dim_1(:, 2) ~= Inf));
    end
    all_Betti_2(i, 1:num_nodes(i)) = struct_4.bettis_orig(3, :);
    bars_dim_2 = struct_4.bars_orig{3, 2};
    if ~isempty(bars_dim_2)
        closed_dim_2(i) = numel(find(bars_dim_2(:, 2) ~= Inf));
    end
    ID(i) = subj_ID;
end

ID = int64(ID); % possibly unnecessary conversion to int type for R

%% Store for analysis in R

clearvars -except ID density num_nodes num_edges ...
    all_C all_d num_change_points all_DoF_C ...
    all_Betti_0 closed_dim_0 ...
    all_Betti_1 closed_dim_1 ...
    all_Betti_2 closed_dim_2

save('/Volumes/My Passport/Curiosity/v8/Data/KNOT/processed_KNOT_data.mat');

