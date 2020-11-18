clc; close all; clear;

%% Compressibility for Growing Knowledge Networks

data_path = '/Volumes/My Passport/Curiosity/Data/Processed';
files = dir(fullfile(data_path, 'subj_*.mat'));

setting = 7;
num_pairs = 100;

% how large is the largest network?
max_time_steps = 0;
for i = 1:length(files)
    load(fullfile(data_path, files(i).name));
    if length(all_adj) > max_time_steps
        max_time_steps = length(all_adj);
    end
end

% initialize empty matrix of compressibility values
C_all = NaN(length(files), max_time_steps);

for i = 1:length(files)
    fprintf('Subject %d of %d.\n', i, length(files))
    load(fullfile(data_path, files(i).name));
    for j = 1:length(all_adj)
        G = full(all_adj{1, j});
        try
            [S, S_low, clusters, Gs] = rate_distortion_upper_info(G, setting, num_pairs);
            C_all(i,j) = mean(S(end) - S);
        catch % assert that invalid return stays NaN
            C_all(i,j) = NaN;
        end
    end
end