clc; close all; clear;

%% Load all necessary data

base_path = '/Volumes/My Passport/Curiosity/';
addpath(genpath(fullfile(base_path, 'Helper')))
addpath(genpath('/Users/sppatankar/Documents/MATLAB/humanStructureFunction'))
data_path = fullfile(base_path, 'v6/Data/KNOT/');
files = dir(fullfile(data_path, 'Raw', 'subj_*.mat'));

% load in 5D scale, C, Betti, DoF
load(fullfile(data_path, 'KNOT_data.mat')); 

density = zeros(1, length(files));
num_nodes = zeros(1, length(files));
num_edges = zeros(1, length(files));


for i = 1:length(files)
    load(fullfile(data_path, 'Raw', files(i).name));
    [density(i), num_nodes(i), num_edges(i)] = density_und(adj);
end

clearvars -except ID DS_5D JE_5D SC_5D ST_5D TS_5D ...
    density num_nodes num_edges ...
    all_C all_DoF all_Betti

save('/Volumes/My Passport/Curiosity/v6/Data/KNOT/KNOT_data_v2.mat');