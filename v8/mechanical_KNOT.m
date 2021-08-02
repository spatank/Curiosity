clc; close all; clear;

%% Preparatory

base_path = '/Volumes/My Passport/Curiosity/';
addpath(genpath(fullfile(base_path, 'Helper')))
addpath(genpath('/Users/sppatankar/Documents/MATLAB/humanStructureFunction'))
data_path = fullfile(base_path, 'v8/Data/KNOT/Preprocessed/');
files = dir(fullfile(data_path, 'subj_*.mat'));

num_iters = 25; % for null models

%% Compile metrics of interest 

for subj_idx = 1:length(files)
        
    fprintf('%s\n', files(subj_idx).name);
    
    load(fullfile(data_path, files(subj_idx).name));
    
    n = size(G, 1);
    d = NaN(1, n);
    DoF = NaN(1, n);
    rigid = NaN(1, n);
    conform = NaN(1, n);
    
    curr_d = 1; 
    
    for i = 1:n
        % fprintf('Nodes %d of %d.\n', i, n);
        G_filt = G(1:i, 1:i);
        % remove nodes of degree zero
        node_degree = sum(G_filt);
        keep_nodes = find(node_degree ~= 0);
        G_filt = G_filt(keep_nodes, keep_nodes);
        if isempty(G_filt)
            continue
        end
        [~, num_nodes, num_edges] = density_und(G_filt);
        curr_DoF = (curr_d * num_nodes) - num_edges; % DoF = rigid + conform
        curr_rigid = curr_d * (curr_d + 1) / 2;
        while curr_DoF - curr_rigid < 0
            curr_d = curr_d + 1;
            curr_DoF = (curr_d * num_nodes) - num_edges; % DoF = rigid + conform
            curr_rigid = curr_d * (curr_d + 1) / 2;
        end
        d(i) = curr_d;
        DoF(i) = curr_DoF;
        rigid(i) = curr_rigid;
        conform(i) = curr_DoF - curr_rigid;
    end

    % Edge rewired networks
    
    d_edge_rewired = NaN(num_iters, n);
    DoF_edge_rewired = NaN(num_iters, n);
    rigid_edge_rewired = NaN(num_iters, n);
    conform_edge_rewired = NaN(num_iters, n);

    for i = 1:num_iters
        
        curr_d = 1;
        
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
                continue
            end
            [~, num_nodes, num_edges] = density_und(G_filt);
            curr_DoF = (curr_d * num_nodes) - num_edges; % DoF = rigid + conform
            curr_rigid = curr_d * (curr_d + 1) / 2;
            while curr_DoF - curr_rigid < 0
                curr_d = curr_d + 1;
                curr_DoF = (curr_d * num_nodes) - num_edges; % DoF = rigid + conform
                curr_rigid = curr_d * (curr_d + 1) / 2;
            end
            d_edge_rewired(i, j) = curr_d;
            DoF_edge_rewired(i, j) = curr_DoF;
            rigid_edge_rewired(i, j) = curr_rigid;
            conform_edge_rewired(i, j) = curr_DoF - curr_rigid;
        end
    end

    % Latticized networks
    
    d_latticized = NaN(num_iters, n);
    DoF_latticized = NaN(num_iters, n);
    rigid_latticized = NaN(num_iters, n);
    conform_latticized = NaN(num_iters, n);

    for i = 1:num_iters
        
        curr_d = 1;
        
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
                continue
            end
            [~, num_nodes, num_edges] = density_und(G_filt);
            curr_DoF = (curr_d * num_nodes) - num_edges; % DoF = rigid + conform
            curr_rigid = curr_d * (curr_d + 1) / 2;
            while curr_DoF - curr_rigid < 0
                curr_d = curr_d + 1;
                curr_DoF = (curr_d * num_nodes) - num_edges; % DoF = rigid + conform
                curr_rigid = curr_d * (curr_d + 1) / 2;
            end
            d_latticized(i, j) = curr_d;
            DoF_latticized(i, j) = curr_DoF;
            rigid_latticized(i, j) = curr_rigid;
            conform_latticized(i, j) = curr_DoF - curr_rigid;
        end
    end

    % Save variables of interest    
    parse_filename = split(files(subj_idx).name, '_');
    subj_ID = str2double(parse_filename{2});
    save([parse_filename{1}, '_', parse_filename{2}, '_mech.mat'], ...
        'd', 'DoF', 'rigid', 'conform', ...
        'd_edge_rewired', 'DoF_edge_rewired', 'rigid_edge_rewired', 'conform_edge_rewired', ...
        'd_latticized', 'DoF_latticized', 'rigid_latticized', 'conform_latticized', ...
        'subj_ID');
end

%% Diagnostic Plotting

% figure;
% hold on
% plot(1:n, d, 'LineWidth', 2, ...
%     'Color', [0, 0, 0]);
% plot(1:n, mean(d_edge_rewired, 'omitnan'), 'LineWidth', 2, ...
%     'Color', [0.8500, 0.3250, 0.0980]);
% plot(1:n, mean(d_latticized, 'omitnan'), 'LineWidth', 2, ...
%     'Color', [0.9290, 0.6940, 0.1250]);
% hold off
% xlabel('Nodes', 'FontSize', 20);
% ylabel('Dimensionality', 'FontSize', 20);
% legend('Original', 'Edges Rewired', 'Latticized', ...
%     'Location', 'SE');
% prettify
% 
% figure;
% hold on
% plot(1:n, rigid, 'LineWidth', 2, ...
%     'Color', [0, 0, 0]);
% plot(1:n, mean(rigid_edge_rewired, 'omitnan'), 'LineWidth', 2, ...
%     'Color', [0.8500, 0.3250, 0.0980]);
% plot(1:n, mean(rigid_latticized, 'omitnan'), 'LineWidth', 2, ...
%     'Color', [0.9290, 0.6940, 0.1250]);
% hold off
% xlabel('Nodes', 'FontSize', 20);
% ylabel('Rigid Body DoF', 'FontSize', 20);
% legend('Original', 'Edges Rewired', 'Latticized', ...
%     'Location', 'SE');
% prettify
% 
% figure;
% hold on
% plot(1:n, conform, 'LineWidth', 2, ...
%     'Color', [0, 0, 0]);
% plot(1:n, mean(conform_edge_rewired, 'omitnan'), 'LineWidth', 2, ...
%     'Color', [0.8500, 0.3250, 0.0980]);
% plot(1:n, mean(conform_latticized, 'omitnan'), 'LineWidth', 2, ...
%     'Color', [0.9290, 0.6940, 0.1250]);
% hold off
% xlabel('Nodes', 'FontSize', 20);
% ylabel('Conformable DoF', 'FontSize', 20);
% legend('Original', 'Edges Rewired', 'Latticized', ...
%     'Location', 'NW');
% prettify
% 
% figure;
% hold on
% plot(1:n, DoF, 'LineWidth', 2, ...
%     'Color', [0, 0, 0]);
% plot(1:n, mean(DoF_edge_rewired, 'omitnan'), 'LineWidth', 2, ...
%     'Color', [0.8500, 0.3250, 0.0980]);
% plot(1:n, mean(DoF_latticized, 'omitnan'), 'LineWidth', 2, ...
%     'Color', [0.9290, 0.6940, 0.1250]);
% hold off
% xlabel('Nodes', 'FontSize', 20);
% ylabel('DoF', 'FontSize', 20);
% legend('Original', 'Edges Rewired', 'Latticized', ...
%     'Location', 'NW');
% prettify

