clc; close all; clear;

%% Preparatory

base_path = '/Volumes/My Passport/Curiosity/v8/Data/KNOT/';
addpath(genpath(fullfile(base_path, 'Helper')))
addpath(genpath('/Users/sppatankar/Documents/MATLAB/humanStructureFunction'))
preprocessed_data_path = fullfile(base_path, 'Preprocessed/');
preprocessed_files = dir(fullfile(preprocessed_data_path, 'subj_*.mat'));
processed_data_path = fullfile(base_path, 'Processed/Mech/');
processed_files = dir(fullfile(processed_data_path, 'subj_*.mat'));

num_iters = 25; % for null models

%% Analysis

% What percentage of the jumps in dimensionality happen when there is no 
% hyperlink between temporally adjacent pages?

% - for each individual
% - IN THE ORIGINAL DATA
% -     collect filtration indices where dimensionality jumps occur
% -     check whether associated traversal was along a hyperlink
% - IN THE NULL MODELS
% -     collect filtration indices where dimensionality jumps occur
% -     check whether associated traversal was along a hyperlink

% What percentage of dimensionality changes in the real data happen
% alongside hyperlinked traversals?
percentage_hyperlinked_changes = zeros(1, length(preprocessed_files));
percentage_hyperlinked_non_changes = zeros(1, length(preprocessed_files));

% What is the percentage in the null data?
percentage_hyperlinked_changes_null = zeros(length(preprocessed_files), num_iters);
percentage_hyperlinked_non_changes_null = zeros(length(preprocessed_files), num_iters);

for i = 1:length(preprocessed_files)

    % load data for adjacency matrices
    load(fullfile(preprocessed_data_path, preprocessed_files(i).name));
    % load data for dimensionality changes
    load(fullfile(processed_data_path, processed_files(i).name));

    orig_change_points = find(diff(d) == 1);
    orig_is_hyperlinked = zeros(size(orig_change_points));
    for j = 1:length(orig_change_points)
        filt_idx = orig_change_points(j);
        if G(filt_idx, filt_idx + 1) == 1
            orig_is_hyperlinked(j) = 1;
        end
    end
    percentage_hyperlinked_changes(i) = ...
        sum(orig_is_hyperlinked)/length(orig_is_hyperlinked);
    
    orig_non_change_points = find(diff(d) ~= 1);
    orig_non_change_points_is_hyperlinked = ...
        zeros(size(orig_non_change_points));
    for k = 1:length(orig_non_change_points)
        filt_idx = orig_non_change_points(k);
        if G(filt_idx, filt_idx + 1) == 1
            orig_non_change_points_is_hyperlinked(k) = 1;
        end
    end
    percentage_hyperlinked_non_changes(i) = ...
        sum(orig_non_change_points_is_hyperlinked)/length(orig_non_change_points_is_hyperlinked);
    
    for iter = 1:num_iters
        G_null = edges_rewired_weighted(:, :, iter);
        G_null(G_null == Inf) = 0;
        G_null(G_null ~= 0) = 1;
        G_null(1:n+1:end) = 0; % remove the diagonal
        null_change_points = find(diff(d_edge_rewired(iter, :)) == 1);
        null_is_hyperlinked = zeros(size(null_change_points));
        for l = 1:length(null_change_points)
            filt_idx = null_change_points(l);
            if G_null(filt_idx, filt_idx + 1) == 1
                null_is_hyperlinked(l) = 1;
            end
        end
        percentage_hyperlinked_changes_null(i, iter) = ...
            sum(null_is_hyperlinked)/length(null_is_hyperlinked);
        
        null_non_change_points = find(diff(d_edge_rewired(iter, :)) ~= 1);
        null_non_change_points_is_hyperlinked = zeros(size(null_non_change_points));
        for l = 1:length(null_non_change_points)
            filt_idx = null_non_change_points(l);
            if G_null(filt_idx, filt_idx + 1) == 1
                null_non_change_points_is_hyperlinked(l) = 1;
            end
        end
        percentage_hyperlinked_non_changes_null(i, iter) = ...
            sum(null_non_change_points_is_hyperlinked)/length(null_non_change_points_is_hyperlinked);
    end
end

mean(percentage_hyperlinked_changes, 'omitnan')


mean(percentage_hyperlinked_non_changes, 'omitnan')


mean(mean(percentage_hyperlinked_changes_null, 'omitnan'), 'omitnan')


mean(mean(percentage_hyperlinked_non_changes_null, 'omitnan'), 'omitnan')