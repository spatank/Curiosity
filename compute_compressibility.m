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

%% Compile Compressibility Values

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

%% Plot Time vs. Compressibility

clc; close all; clear;

load('/Volumes/My Passport/Curiosity/Data/C_all.mat')

% some compressibility values are complex
C_all(imag(C_all) ~= 0) = NaN;
% some are implausibly large positive/negative numbers
C_all(abs(C_all) > 100) = NaN;

time_steps = 1:length(C_all);
C_mean = mean(C_all, 'omitnan');

f = figure('color', 'w');
scatter(time_steps, C_mean, 'filled');
xlabel('Time Steps', 'FontSize', 15);
ylabel('Compressibility', 'FontSize', 15);
title('Time vs. Compressibility', 'FontSize', 15);


