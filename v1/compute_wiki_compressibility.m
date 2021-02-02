clc; clear; close all;

%% Compressibility for Growing Knowledge Networks

topic = 'physics';
base_data_path = '/Volumes/My Passport/Curiosity/Data/wiki_data_processed/';
full_data_path = fullfile(base_data_path, topic);
files = dir(fullfile(full_data_path, strcat(topic, '*.mat')));

timeline = zeros(1, length(files));
compressibility = zeros(1, length(files));
num_components = zeros(1, length(files));

setting = 7;
num_pairs = 100;

for i = 1:length(files)
    fprintf('Year %d of %d.\n', i, length(files));
    load(fullfile(full_data_path, files(i).name));
    timeline(i) = year;
    G = full(adj);
    % G(G > 0) = 1; % binarize
    G = (G + G') - eye(size(G, 1)) .* diag(G); % make symmetric
    if isempty(G)
        fprintf('Year %d of %d is empty.\n', i, length(files));
        compressibility(i) = NaN;
        continue
    end
    [components, component_sizes] = conncomp(digraph(G), 'Type', 'Weak');
    num_components(i) = max(components);
    idx = component_sizes(components) == max(component_sizes);
    largest_G = full(adjacency(subgraph(digraph(G), idx)));
    try
        [S, S_low, clusters, Gs] = rate_distortion_upper_info_new(largest_G, setting, num_pairs);
        compressibility(i) = mean(S(end) - S);
    catch
        compressibility(i) = NaN;
    end
end

%% Plot

f = figure('color', 'w');
scatter(timeline, compressibility, 'b', 'filled');
xlabel('Year', 'FontSize', 15);
ylabel('Compressibility', 'FontSize', 15);
% title('Evolutionary Biology Growing Graph', 'FontSize', 15);

f = figure('color', 'w');
scatter(timeline, num_components, 'r', 'filled');
xlabel('Year', 'FontSize', 15);
ylabel('Connected Components', 'FontSize', 15);
% title('Evolutionary Biology Growing Graph', 'FontSize', 15);

% for i = 1:length(edge_info)
%     fprintf('from: %s, to: %s.\n', edge_info{1, i}.from, edge_info{1, i}.to);
% end