clc; close all; clear;

%% Preparatory

base_path = '/Volumes/My Passport/Curiosity/';
addpath(genpath(fullfile(base_path, 'Helper')))
addpath(genpath('/Users/sppatankar/Documents/MATLAB/humanStructureFunction'))
data_path = fullfile(base_path, 'v6/Data/KNOT/Preprocessed/');
files = dir(fullfile(data_path, 'subj_*.mat'));

num_iters = 25; % for null models

%% Compile metrics of interest 

% for subj_idx = 20:length(files)
for subj_idx = 1:length(files)
        
    fprintf('%s\n', files(subj_idx).name);
    
    load(fullfile(data_path, files(subj_idx).name));
    
    n = size(G, 1);
    d = zeros(1, n);
    
    for i = 1:n
        % fprintf('Nodes %d of %d.\n', i, n);
        G_filt = G(1:i, 1:i);
        % remove nodes of degree zero
        node_degree = sum(G_filt);
        keep_nodes = find(node_degree ~= 0);
        G_filt = G_filt(keep_nodes, keep_nodes);
        if isempty(G_filt)
            d(i) = NaN;
            continue
        end
        [~, num_nodes, num_edges] = density_und(G_filt);
        p = [1, (1 - (2*num_nodes)), 2*num_edges];
        r = roots(p);
        r = r(r > 0); % remove negative roots
%         r = ceil(min(r));
        r = (min(r));
        d(i) = r;
    end

    % Edge rewired networks

    d_edge_rewired = zeros(num_iters, n);

    for i = 1:num_iters
        curr_weighted_G = edges_rewired_weighted(:, :, i);
        curr_weighted_G(curr_weighted_G == Inf) = 0;
        curr_weighted_G(curr_weighted_G ~= 0) = 1;
        curr_weighted_G(1:n+1:end) = 0; % remove the diagonal
        curr_G = curr_weighted_G;
        symmetry_flag = issymmetric(curr_G);
        % fprintf('Iteration %d of %d; symmetry: %d.\n', i, num_iters, symmetry_flag);
        for j = 1:n
            G_filt = curr_G(1:j, 1:j);
            % remove nodes of degree zero
            node_degree = sum(G_filt);
            keep_nodes = find(node_degree ~= 0);
            G_filt = G_filt(keep_nodes, keep_nodes);
            if isempty(G_filt)
                d_edge_rewired(i, j) = NaN;
                continue
            end
            [~, num_nodes, num_edges] = density_und(G_filt);
            p = [1, (1 - (2*num_nodes)), 2*num_edges];
            r = roots(p);
            r = r(r > 0); % remove negative roots
%             r = ceil(min(r));
            r = (min(r));
            d_edge_rewired(i, j) = r;
        end
    end

    % Latticized networks

    d_latticized = zeros(num_iters, n);

    for i = 1:num_iters
        curr_weighted_G = latticized_weighted(:, :, i);
        curr_weighted_G(curr_weighted_G == Inf) = 0;
        curr_weighted_G(curr_weighted_G ~= 0) = 1;
        curr_weighted_G(1:n+1:end) = 0;
        curr_G = curr_weighted_G;
        symmetry_flag = issymmetric(curr_G);
        % fprintf('Iteration %d of %d; symmetry: %d.\n', i, num_iters, symmetry_flag);
        for j = 1:n
            G_filt = curr_G(1:j, 1:j);
            % remove nodes of degree zero
            node_degree = sum(G_filt);
            keep_nodes = find(node_degree ~= 0);
            G_filt = G_filt(keep_nodes, keep_nodes);
            if isempty(G_filt)
                d_latticized(i, j) = NaN;
                continue
            end
            [~, num_nodes, num_edges] = density_und(G_filt);
            p = [1, (1 - (2*num_nodes)), 2*num_edges];
            r = roots(p);
            r = r(r > 0); % remove negative roots
%             r = ceil(min(r));
            r = (min(r));
            d_latticized(i, j) = r;
        end
    end

    % Save variables of interest    
    parse_filename = split(files(subj_idx).name, '_');
    subj_ID = str2double(parse_filename{2});
    save([parse_filename{1}, '_', parse_filename{2}, '_d.mat'], ...
        'd', 'd_edge_rewired', 'd_latticized', 'n', 'subj_ID');
end