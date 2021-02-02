clc; close all; clear;

%% Shared Parameters

setting = 7;
num_pairs = 100;
base_data_path = '/Volumes/My Passport/Curiosity/Data/SimulatedNetworks';

%% Probabilistic Triangle Model

ps_PT = 0.3:0.2:0.9;
thresholds = 0.1:0.1:1.0;
all_compressibilities_PT = zeros(length(ps_PT), length(thresholds));

for i = 1:length(ps_PT)
    p = ps_PT(i);
    files = dir(fullfile(base_data_path, 'Triangle', num2str(p), 'p*.mat'));
    curr_p_compressibilities = zeros(length(files), length(thresholds));
    for j = 1:length(files)
        fprintf('PT, p = %0.2f, iter = %d\n', p, j);
        load(fullfile(base_data_path, 'Triangle', num2str(p), files(j).name))
        % load in iter, p, thresholded_Gs
        for k = 1:length(thresholds)
            G = thresholded_Gs(:, :, k);
            [components, component_sizes] = conncomp(digraph(G), 'Type', 'Weak');
            idx = component_sizes(components) == max(component_sizes);
            largest_G = full(adjacency(subgraph(digraph(G), idx)));
            try
                [S, S_low, clusters, Gs] = rate_distortion_upper_info_new(largest_G, setting, num_pairs);
                curr_p_compressibilities(j, k) = mean(S(end) - S);
            catch
                curr_p_compressibilities(j, k) = NaN;
            end
        end
    end
    all_compressibilities_PT(i, :) = mean(curr_p_compressibilities, 'omitnan');
end


%% Weighted Probabilistic Triangle Model

ps = 0.1:0.2:0.9;
thresholds = 0.1:0.1:1.0;
all_compressibilities_WPT = zeros(length(ps), length(thresholds));

for i = 1:length(ps)
    p = ps(i);
    files = dir(fullfile(base_data_path, 'WeightedTriangle', num2str(p), 'p*.mat'));
    curr_p_compressibilities = zeros(length(files), length(thresholds));
    for j = 1:length(files)
        fprintf('WPT, p = %0.2f, iter = %d\n', p, j);
        load(fullfile(base_data_path, 'WeightedTriangle', num2str(p), files(j).name))
        % load in iter, p, thresholded_Gs
        for k = 1:length(thresholds)
            G = thresholded_Gs(:, :, k);
            [components, component_sizes] = conncomp(digraph(G), 'Type', 'Weak');
            idx = component_sizes(components) == max(component_sizes);
            largest_G = full(adjacency(subgraph(digraph(G), idx)));
            try
                [S, S_low, clusters, Gs] = rate_distortion_upper_info_new(largest_G, setting, num_pairs);
                curr_p_compressibilities(j, k) = mean(S(end) - S);
            catch
                curr_p_compressibilities(j, k) = NaN;
            end
        end
    end
    all_compressibilities_WPT(i, :) = mean(curr_p_compressibilities, 'omitnan');
end


%% IID

clearvars -except setting num_pairs base_data_path ps ps_PT thresholds all_compressibilities_PT all_compressibilities_WPT 

files = dir(fullfile(base_data_path, 'IID', 'iid*.mat'));
all_compressibilies_IID = zeros(length(files), length(thresholds));
for i = 1:length(files)
    fprintf('IID iter = %d\n', i);
    load(fullfile(base_data_path, 'IID', files(i).name))
    % load in iter, p, thresholded_Gs
    for j = 1:length(thresholds)
        G = thresholded_Gs(:, :, j);
        [components, component_sizes] = conncomp(digraph(G), 'Type', 'Weak');
        idx = component_sizes(components) == max(component_sizes);
        largest_G = full(adjacency(subgraph(digraph(G), idx)));
        try
            [S, S_low, clusters, Gs] = rate_distortion_upper_info_new(largest_G, setting, num_pairs);
            all_compressibilies_IID(i, j) = mean(S(end) - S);
        catch
            all_compressibilies_IID(i, j) = NaN;
        end
    end
end
all_compressibilities_IID = mean(all_compressibilies_IID, 'omitnan');


%% Plot

close all;

ps_PT = 0.3:0.2:0.9;

figure
hold on
for i = 1:length(ps_PT)
    plot(thresholds, all_compressibilities_PT(i, :), 'LineWidth', 2);
end
plot(thresholds, all_compressibilities_IID, 'LineWidth', 2, 'LineStyle', '--');
xlabel('\rho', 'FontSize', 20);
ylabel('Compressibility', 'FontSize', 15);
title('Probabilistic Triangle', 'FontSize', 15)
legend('0.3', '0.5', '0.7', '0.9', 'IID', 'Location', 'SouthEast')
hold off
prettify

figure
hold on
for i = 1:length(ps)
    plot(thresholds, all_compressibilities_WPT(i, :), 'LineWidth', 2);
end
xlabel('\rho', 'FontSize', 20);
ylabel('Compressibility', 'FontSize', 15);
title('Weighted Probabilistic Triangle', 'FontSize', 15)
legend('0.1', '0.3', '0.5', '0.7', '0.9', 'Location', 'SouthEast')
hold off
prettify

figure
plot(thresholds, all_compressibilities_IID, 'LineWidth', 2);
xlabel('\rho', 'FontSize', 20);
ylabel('Compressibility', 'FontSize', 15);
title('IID', 'FontSize', 15);
prettify

figure
hold on
plot(thresholds, all_compressibilities_WPT(1, :), 'LineWidth', 2);
plot(thresholds, all_compressibilities_WPT(3, :), 'LineWidth', 2);
plot(thresholds, all_compressibilities_WPT(5, :), 'LineWidth', 2);
plot(thresholds, all_compressibilities_IID, 'LineWidth', 2, 'LineStyle', '--');
xlabel('\rho', 'FontSize', 20);
ylabel('Compressibility', 'FontSize', 15);
legend('0.1', '0.5', '0.9', 'IID', 'Location', 'SouthEast')
hold off
prettify

