clc; close all; clear;

%% Preparatory

base_path = '/Volumes/My Passport/Curiosity/';
addpath(genpath(fullfile(base_path, 'Helper')))
addpath(genpath('/Users/sppatankar/Documents/MATLAB/humanStructureFunction'))
data_path_C_DoF = fullfile(base_path, 'v6/Data/KNOT/Processed/C_DoF/');
files_C_DoF = dir(fullfile(data_path_C_DoF, 'subj_*.mat'));
data_path_Betti = fullfile(base_path, 'v6/Data/KNOT/Processed/Betti/');
files_Betti = dir(fullfile(data_path_Betti, 'subj_*.mat'));

% how large is the largest network?
max_size = 0;
for i = 1:length(files_Betti)
    load(fullfile(data_path_Betti, files_Betti(i).name));
    if length(bettis_orig) > max_size
        max_size = length(bettis_orig);
    end
end

%% Collect and store metrics

all_C = NaN(length(files_C_DoF), max_size);
all_DoF = NaN(length(files_C_DoF), max_size);
all_Betti = NaN(length(files_C_DoF), max_size);

for i = 1:length(files_C_DoF)
    load(strcat(data_path_C_DoF, files_C_DoF(i).name));
    load(strcat(data_path_Betti, files_Betti(i).name));
    all_C(i, 1:length(C)) = C;
    all_DoF(i, 1:length(DoF)) = DoF;
    all_Betti(i, 1:length(bettis_orig)) = bettis_orig(2, :);
end

%% Store for analysis in R

clearvars -except all_C all_DoF all_Betti
load('/Volumes/My Passport/Curiosity/Data/five_D.mat')
save('/Volumes/My Passport/Curiosity/v6/Data/KNOT/KNOT_data.mat');