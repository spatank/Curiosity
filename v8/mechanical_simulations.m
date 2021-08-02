clc; close all; clear;

%% Preparatory

base_path = '/Volumes/My Passport/Curiosity/';
addpath(genpath(fullfile(base_path, 'Helper')))
data_path = fullfile(base_path, 'v8/Data/Simulations/');
load(fullfile(data_path, 'simulated_nets.mat'));

num_iters = 100; % for null models
n = 70; % number of nodes in simulated networks

%% Compile metrics of interest 

d_CP = NaN(num_iters, n);
DoF_CP = NaN(num_iters, n);
rigid_CP = NaN(num_iters, n);
conform_CP = NaN(num_iters, n);

d_PP = NaN(num_iters, n);
DoF_PP = NaN(num_iters, n);
rigid_PP = NaN(num_iters, n);
conform_PP = NaN(num_iters, n);

d_PA = NaN(num_iters, n);
DoF_PA = NaN(num_iters, n);
rigid_PA = NaN(num_iters, n);
conform_PA = NaN(num_iters, n);

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
        d_CP(iter, i) = curr_d;
        DoF_CP(iter, i) = curr_DoF;
        rigid_CP(iter, i) = curr_rigid;
        conform_CP(iter, i) = curr_DoF - curr_rigid;
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
        d_PP(iter, i) = curr_d;
        DoF_PP(iter, i) = curr_DoF;
        rigid_PP(iter, i) = curr_rigid;
        conform_PP(iter, i) = curr_DoF - curr_rigid;
    end
    
    % Preferential attachment networks
    G = PA_nets(:, :, iter);
    
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
        d_PA(iter, i) = curr_d;
        DoF_PA(iter, i) = curr_DoF;
        rigid_PA(iter, i) = curr_rigid;
        conform_PA(iter, i) = curr_DoF - curr_rigid;
    end
end

clearvars -except n ...
    d_CP DoF_CP rigid_CP conform_CP ...
    d_PP DoF_PP rigid_PP conform_PP ...
    d_PA DoF_PA rigid_PA conform_PA

% Save variables of interest    
filename = 'simulated_nets_mech.mat';
save(filename, 'd_CP', 'DoF_CP', 'rigid_CP', 'conform_CP', ...
    'd_PP', 'DoF_PP', 'rigid_PP', 'conform_PP', ...
    'd_PA', 'DoF_PA', 'rigid_PA', 'conform_PA');

%% Diagnostic Plotting

% Constant probability model

figure;
subplot(2, 2, 1)
hold on
plot(1:n, d_CP, 'LineWidth', 2, ...
    'Color', [0.7, 0.7, 0.7]);
plot(1:n, mean(d_CP, 'omitnan'), 'LineWidth', 2, ...
    'Color', [0, 0, 0]);
hold off
xlabel('Nodes', 'FontSize', 20);
ylabel('Dimensionality', 'FontSize', 20);
title('CP', 'FontSize', 20);
prettify

% figure;
subplot(2, 2, 2)
hold on
plot(1:n, rigid_CP, 'LineWidth', 2, ...
    'Color', [0.7, 0.7, 0.7]);
plot(1:n, mean(rigid_CP, 'omitnan'), 'LineWidth', 2, ...
    'Color', [0, 0, 0]);
hold off
xlabel('Nodes', 'FontSize', 20);
ylabel('Rigid Body DoF', 'FontSize', 20);
title('CP', 'FontSize', 20);
prettify

% figure;
subplot(2, 2, 3)
hold on
plot(1:n, conform_CP, 'LineWidth', 2, ...
    'Color', [0.7, 0.7, 0.7]);
plot(1:n, mean(conform_CP, 'omitnan'), 'LineWidth', 2, ...
    'Color', [0, 0, 0]);
hold off
xlabel('Nodes', 'FontSize', 20);
ylabel('Conformable DoF', 'FontSize', 20);
title('CP', 'FontSize', 20);
prettify

% figure;
subplot(2, 2, 4)
hold on
plot(1:n, DoF_CP, 'LineWidth', 2, ...
    'Color', [0.7, 0.7, 0.7]);
plot(1:n, mean(DoF_CP, 'omitnan'), 'LineWidth', 2, ...
    'Color', [0, 0, 0]);
hold off
xlabel('Nodes', 'FontSize', 20);
ylabel('DoF', 'FontSize', 20);
title('CP', 'FontSize', 20);
prettify

% Proportional probability

figure;
subplot(2, 2, 1)
hold on
plot(1:n, d_PP, 'LineWidth', 2, ...
    'Color', [0.7, 0.7, 0.7]);
plot(1:n, mean(d_PP, 'omitnan'), 'LineWidth', 2, ...
    'Color', [0, 0, 0]);
hold off
xlabel('Nodes', 'FontSize', 20);
ylabel('Dimensionality', 'FontSize', 20);
title('PP', 'FontSize', 20);
prettify

% figure;
subplot(2, 2, 2)
hold on
plot(1:n, rigid_PP, 'LineWidth', 2, ...
    'Color', [0.7, 0.7, 0.7]);
plot(1:n, mean(rigid_PP, 'omitnan'), 'LineWidth', 2, ...
    'Color', [0, 0, 0]);
hold off
xlabel('Nodes', 'FontSize', 20);
ylabel('Rigid Body DoF', 'FontSize', 20);
title('PP', 'FontSize', 20);
prettify

% figure;
subplot(2, 2, 3)
hold on
plot(1:n, conform_PP, 'LineWidth', 2, ...
    'Color', [0.7, 0.7, 0.7]);
plot(1:n, mean(conform_PP, 'omitnan'), 'LineWidth', 2, ...
    'Color', [0, 0, 0]);
hold off
xlabel('Nodes', 'FontSize', 20);
ylabel('Conformable DoF', 'FontSize', 20);
title('PP', 'FontSize', 20);
prettify

% figure;
subplot(2, 2, 4)
hold on
plot(1:n, DoF_PP, 'LineWidth', 2, ...
    'Color', [0.7, 0.7, 0.7]);
plot(1:n, mean(DoF_PP, 'omitnan'), 'LineWidth', 2, ...
    'Color', [0, 0, 0]);
hold off
xlabel('Nodes', 'FontSize', 20);
ylabel('DoF', 'FontSize', 20);
title('PP', 'FontSize', 20);
prettify

% Preferential attachment model

figure;
subplot(2, 2, 1)
hold on
plot(1:n, d_PA, 'LineWidth', 2, ...
    'Color', [0.7, 0.7, 0.7]);
plot(1:n, mean(d_PA, 'omitnan'), 'LineWidth', 2, ...
    'Color', [0, 0, 0]);
hold off
xlabel('Nodes', 'FontSize', 20);
ylabel('Dimensionality', 'FontSize', 20);
title('PA', 'FontSize', 20);
prettify

% figure;
subplot(2, 2, 2)
hold on
plot(1:n, rigid_PA, 'LineWidth', 2, ...
    'Color', [0.7, 0.7, 0.7]);
plot(1:n, mean(rigid_PA, 'omitnan'), 'LineWidth', 2, ...
    'Color', [0, 0, 0]);
hold off
xlabel('Nodes', 'FontSize', 20);
ylabel('Rigid Body DoF', 'FontSize', 20);
title('PA', 'FontSize', 20);
prettify

% figure;
subplot(2, 2, 3)
hold on
plot(1:n, conform_PA, 'LineWidth', 2, ...
    'Color', [0.7, 0.7, 0.7]);
plot(1:n, mean(conform_PA, 'omitnan'), 'LineWidth', 2, ...
    'Color', [0, 0, 0]);
hold off
xlabel('Nodes', 'FontSize', 20);
ylabel('Conformable DoF', 'FontSize', 20);
title('PA', 'FontSize', 20);
prettify

% figure;
subplot(2, 2, 4)
hold on
plot(1:n, DoF_PA, 'LineWidth', 2, ...
    'Color', [0.7, 0.7, 0.7]);
plot(1:n, mean(DoF_PA, 'omitnan'), 'LineWidth', 2, ...
    'Color', [0, 0, 0]);
hold off
xlabel('Nodes', 'FontSize', 20);
ylabel('DoF', 'FontSize', 20);
title('PA', 'FontSize', 20);
prettify

