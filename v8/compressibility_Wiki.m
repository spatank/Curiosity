clc; close all; clear;

%% Preparatory

base_path = '/Volumes/My Passport/Curiosity/';
addpath(genpath(fullfile(base_path, 'Helper')))
addpath(genpath('/Users/sppatankar/Documents/MATLAB/humanStructureFunction'))
data_path = fullfile(base_path, 'v7/Data/Wiki/Preprocessed/');
files = dir(fullfile(data_path, '*_preprocessed.mat'));

setting = 7;
num_pairs = 100;
num_iters = 25; % for null models

%% Compile metrics of interest 

for topic_idx = 1:length(files)
        
    fprintf('%s\n', files(topic_idx).name);
    
    load(fullfile(data_path, files(topic_idx).name));
    
    n = size(G, 1);
    clust_coef = NaN(1, n);
    C = NaN(1, n);
    C_norm = NaN(1, n);
    entropy = NaN(1, n);
    
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
            C_norm(i) = C(i)/S(end);
            entropy(i) = S(end);
        catch
            fprintf('Error caught; size. %d\n', i);
        end
    end

    % Edge rewired networks
    
    clust_coef_edge_rewired = NaN(num_iters, n);
    C_edge_rewired = NaN(num_iters, n);
    C_norm_edge_rewired = NaN(num_iters, n);
    entropy_edge_rewired = NaN(num_iters, n);
    

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
                C_norm_edge_rewired(i, j) = C_edge_rewired(i, j)/S(end);
                entropy_edge_rewired(i, j) = S(end);
            catch
                fprintf('Rewired: error in iter %d at stage %d.\n', i, j);
            end
        end
    end

    % Latticized networks
    
    clust_coef_latticized = NaN(num_iters, n);
    C_latticized = NaN(num_iters, n);
    C_norm_latticized = NaN(num_iters, n);
    entropy_latticized = NaN(num_iters, n);
    

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
                C_norm_latticized(i, j) = C_latticized(i, j)/S(end);
                entropy_latticized(i, j) = S(end);
            catch
                fprintf('Latticized: error in iter %d at stage %d.\n', i, j);
            end
        end
    end

    % Save variables of interest    
    parse_filename = split(files(topic_idx).name, '_');
    topic_ID = strjoin(parse_filename(1:end-1), '_');
    save(strcat(topic_ID, '_C.mat'), ...
        'clust_coef', 'C', 'C_norm', 'entropy', ...
        'clust_coef_edge_rewired', 'C_edge_rewired', 'C_norm_edge_rewired', 'entropy_edge_rewired', ...
        'clust_coef_latticized', 'C_latticized', 'C_norm_latticized', 'entropy_latticized', ...
        'n', 'topic_ID');
    
end

