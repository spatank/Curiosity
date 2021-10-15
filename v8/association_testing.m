clc; close all; clear;

%% Preparatory

base_path = '/Volumes/My Passport/Curiosity/';

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

ID = NaN(length(files_PH), 1);
filled_dim_0_cycles = NaN(length(files_PH), 1);
filled_dim_1_cycles = NaN(length(files_PH), 1);
filled_dim_2_cycles = NaN(length(files_PH), 1);

for i = 1:length(files_PH)
    filename = files_PH(i).name;
    subj_ID = strsplit(filename, '_');
    subj_ID = str2double(cell2mat(subj_ID(2))); % 1: `subj', 2: ID, 3: 'PH.mat'
    ID(i) = subj_ID;
    data_struct = load(strcat(data_path_PH, filename));
    num_nodes = double(num_nodes); % int64 to double for downstream division
    bars_dim_0 = data_struct.bars_orig{1, 2};
    if ~isempty(bars_dim_0)
        filled_dim_0_cycles(i) = numel(find(bars_dim_0(:, 2) ~= Inf));
    end
    bars_dim_1 = data_struct.bars_orig{2, 2};
    if ~isempty(bars_dim_1)
        filled_dim_1_cycles(i) = numel(find(bars_dim_1(:, 2) ~= Inf));
    end
    bars_dim_2 = data_struct.bars_orig{3, 2};
    if ~isempty(bars_dim_2)
        filled_dim_2_cycles(i) = numel(find(bars_dim_2(:, 2) ~= Inf));
    end
end

new_data = table(ID, ...
    filled_dim_0_cycles, filled_dim_1_cycles, filled_dim_2_cycles);
old_data = readtable('/Volumes/My Passport/Curiosity/v8/Data/KNOT/processed_KNOT_data_no_curves.csv');
all_data = join(old_data, new_data);

clearvars -except all_data