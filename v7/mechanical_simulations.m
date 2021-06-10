clc; close all; clear;

%% Preparatory

base_path = '/Volumes/My Passport/Curiosity/';
addpath(genpath(fullfile(base_path, 'Helper')))
data_path = fullfile(base_path, 'v7/Data/Simulations/Preprocessed/');
load(fullfile(data_path, 'simulated_nets.mat'));

num_iters = 50; % for null models
n = 75; % number of nodes in simulated networks

%% Compile metrics of interest 

d_const = NaN(num_iters, n);
DoF_const = NaN(num_iters, n);
rigid_const = NaN(num_iters, n);
conform_const = NaN(num_iters, n);

d_prop = NaN(num_iters, n);
DoF_prop = NaN(num_iters, n);
rigid_prop = NaN(num_iters, n);
conform_prop = NaN(num_iters, n);

for iter = 1:num_iters
        
    fprintf('Iteration %d.\n', iter);
    
    % Constant probability networks
    G = const_prob_nets(:, :, iter);
    
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
        d_const(iter, i) = curr_d;
        DoF_const(iter, i) = curr_DoF;
        rigid_const(iter, i) = curr_rigid;
        conform_const(iter, i) = curr_DoF - curr_rigid;
    end
    
    % Proportional probability networks
    G = prop_prob_nets(:, :, iter);
    
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
        d_prop(iter, i) = curr_d;
        DoF_prop(iter, i) = curr_DoF;
        rigid_prop(iter, i) = curr_rigid;
        conform_prop(iter, i) = curr_DoF - curr_rigid;
    end
end

clearvars -except n d_const DoF_const rigid_const conform_const ...
    d_prop DoF_prop rigid_prop conform_prop

% Save variables of interest    
filename = 'simulated_nets_mech.mat';
save(filename, 'd_const', 'DoF_const', 'rigid_const', 'conform_const', ...
    'd_prop', 'DoF_prop', 'rigid_prop', 'conform_prop');

%% Diagnostic Plotting

figure;
hold on
plot(1:n, d_const, 'LineWidth', 2, ...
    'Color', [0.7, 0.7, 0.7]);
plot(1:n, mean(d_const, 'omitnan'), 'LineWidth', 2, ...
    'Color', [0, 0, 0]);
hold off
xlabel('Nodes', 'FontSize', 20);
ylabel('Dimensionality', 'FontSize', 20);
title('CP', 'FontSize', 20);
prettify

figure;
hold on
plot(1:n, rigid_const, 'LineWidth', 2, ...
    'Color', [0.7, 0.7, 0.7]);
plot(1:n, mean(rigid_const, 'omitnan'), 'LineWidth', 2, ...
    'Color', [0, 0, 0]);
hold off
xlabel('Nodes', 'FontSize', 20);
ylabel('Rigid Body DoF', 'FontSize', 20);
title('CP', 'FontSize', 20);
prettify

figure;
hold on
plot(1:n, conform_const, 'LineWidth', 2, ...
    'Color', [0.7, 0.7, 0.7]);
plot(1:n, mean(conform_const, 'omitnan'), 'LineWidth', 2, ...
    'Color', [0, 0, 0]);
hold off
xlabel('Nodes', 'FontSize', 20);
ylabel('Conformable DoF', 'FontSize', 20);
title('CP', 'FontSize', 20);
prettify

figure;
hold on
plot(1:n, DoF_const, 'LineWidth', 2, ...
    'Color', [0.7, 0.7, 0.7]);
plot(1:n, mean(DoF_const, 'omitnan'), 'LineWidth', 2, ...
    'Color', [0, 0, 0]);
hold off
xlabel('Nodes', 'FontSize', 20);
ylabel('DoF', 'FontSize', 20);
title('CP', 'FontSize', 20);
prettify

figure;
hold on
plot(1:n, d_prop, 'LineWidth', 2, ...
    'Color', [0.7, 0.7, 0.7]);
plot(1:n, mean(d_prop, 'omitnan'), 'LineWidth', 2, ...
    'Color', [0, 0, 0]);
hold off
xlabel('Nodes', 'FontSize', 20);
ylabel('Dimensionality', 'FontSize', 20);
title('PP', 'FontSize', 20);
prettify

figure;
hold on
plot(1:n, rigid_prop, 'LineWidth', 2, ...
    'Color', [0.7, 0.7, 0.7]);
plot(1:n, mean(rigid_prop, 'omitnan'), 'LineWidth', 2, ...
    'Color', [0, 0, 0]);
hold off
xlabel('Nodes', 'FontSize', 20);
ylabel('Rigid Body DoF', 'FontSize', 20);
title('PP', 'FontSize', 20);
prettify

figure;
hold on
plot(1:n, conform_prop, 'LineWidth', 2, ...
    'Color', [0.7, 0.7, 0.7]);
plot(1:n, mean(conform_prop, 'omitnan'), 'LineWidth', 2, ...
    'Color', [0, 0, 0]);
hold off
xlabel('Nodes', 'FontSize', 20);
ylabel('Conformable DoF', 'FontSize', 20);
title('PP', 'FontSize', 20);
prettify

figure;
hold on
plot(1:n, DoF_prop, 'LineWidth', 2, ...
    'Color', [0.7, 0.7, 0.7]);
plot(1:n, mean(DoF_prop, 'omitnan'), 'LineWidth', 2, ...
    'Color', [0, 0, 0]);
hold off
xlabel('Nodes', 'FontSize', 20);
ylabel('DoF', 'FontSize', 20);
title('PP', 'FontSize', 20);
prettify

