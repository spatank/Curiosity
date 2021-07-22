clc; close all; clear;

%% Preparatory

base_path = '/Volumes/My Passport/Curiosity/';
addpath(genpath(fullfile(base_path, 'Helper')))
addpath(genpath('/Users/sppatankar/Documents/MATLAB/humanStructureFunction'))
data_path = fullfile(base_path, 'v7/Data/KNOT/Preprocessed/');
files = dir(fullfile(data_path, 'subj_*.mat'));

setting = 7;
num_pairs = 100;
num_iters = 25; % for null models

%% Compile metrics of interest 

for subj_idx = 1:length(files)
        
    fprintf('%s\n', files(subj_idx).name);
    
    load(fullfile(data_path, files(subj_idx).name));
    
    n = size(G, 1);
    clust_coef = NaN(1, n);
    C = NaN(1, n);
    
    for i = 1:n
        fprintf('Nodes %d of %d.\n', i, n);
        G_filt = G(1:i, 1:i);
        % remove nodes of degree zero
        node_degree = sum(G_filt);
        keep_nodes = find(node_degree ~= 0);
        G_filt = G_filt(keep_nodes, keep_nodes);
        if isempty(G_filt)
            continue
        end
        clust_coef(i) = mean(clustering_coef_bu(G_filt));
        try
            [S, S_low, clusters, Gs] = ...
                rate_distortion_upper_info(G_filt, setting, num_pairs);
            C(i) = mean(S(end) - S);
        catch
            fprintf('Error caught; size. %d\n', i);
            C(i) = NaN;
        end
    end

    % Edge rewired networks

    clust_coef_edge_rewired = NaN(num_iters, n);
    C_edge_rewired = NaN(num_iters, n);
    
    for i = 1:num_iters
        curr_weighted_G = edges_rewired_weighted(:, :, i);
        curr_weighted_G(curr_weighted_G == Inf) = 0;
        curr_weighted_G(curr_weighted_G ~= 0) = 1;
        curr_weighted_G(1:n+1:end) = 0; % remove the diagonal
        curr_G = curr_weighted_G;
        symmetry_flag = issymmetric(curr_G);
        fprintf('Iteration %d of %d; symmetry: %d.\n', i, num_iters, symmetry_flag);
        for j = 1:n
            G_filt = curr_G(1:j, 1:j);
            % remove nodes of degree zero
            node_degree = sum(G_filt);
            keep_nodes = find(node_degree ~= 0);
            G_filt = G_filt(keep_nodes, keep_nodes);
            if isempty(G_filt)
                continue
            end
            clust_coef_edge_rewired(i, j) = mean(clustering_coef_bu(G_filt));
            try
                [S, S_low, clusters, Gs] = ...
                    rate_distortion_upper_info(G_filt, setting, num_pairs);
                C_edge_rewired(i, j) = mean(S(end) - S);
            catch
                C_edge_rewired(i, j) = NaN;
                fprintf('Rewired: error in iter %d at stage %d.\n', i, j);
            end
        end
    end

    % Latticized networks

    clust_coef_latticized = NaN(num_iters, n);
    C_latticized = NaN(num_iters, n);
    
    for i = 1:num_iters
        curr_weighted_G = latticized_weighted(:, :, i);
        curr_weighted_G(curr_weighted_G == Inf) = 0;
        curr_weighted_G(curr_weighted_G ~= 0) = 1;
        curr_weighted_G(1:n+1:end) = 0;
        curr_G = curr_weighted_G;
        symmetry_flag = issymmetric(curr_G);
        fprintf('Iteration %d of %d; symmetry: %d.\n', i, num_iters, symmetry_flag);
        for j = 1:n
            G_filt = curr_G(1:j, 1:j);
            % remove nodes of degree zero
            node_degree = sum(G_filt);
            keep_nodes = find(node_degree ~= 0);
            G_filt = G_filt(keep_nodes, keep_nodes);
            if isempty(G_filt)
                continue
            end
            clust_coef_latticized(i, j) = mean(clustering_coef_bu(G_filt));
            try
                [S, S_low, clusters, Gs] = ...
                    rate_distortion_upper_info(G_filt, setting, num_pairs);
                C_latticized(i, j) = mean(S(end) - S);
            catch
                C_latticized(i, j) = NaN;
                fprintf('Latticized: error in iter %d at stage %d.\n', i, j);
            end
        end
    end

    % Save variables of interest    
    parse_filename = split(files(subj_idx).name, '_');
    subj_ID = str2double(parse_filename{2});
    save([parse_filename{1}, '_', parse_filename{2}, '_C_clust_coef.mat'], ...
        'clust_coef', 'C', 'clust_coef_edge_rewired', 'clust_coef_latticized', ...
        'C_edge_rewired', 'C_latticized', 'n', 'subj_ID');
    
end