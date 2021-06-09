clc; close all; clear;

%% Specify paths and miscellaneous settings

base_path = '/Volumes/My Passport/Curiosity/';
addpath(genpath(fullfile(base_path, 'Helper')))

%% Build networks

n = 75; % number of nodes
p = 0.5; % parameter for constant probability model
iters = 50;

const_prob_nets = NaN(n, n, iters);
const_prob_nets_weighted = NaN(n, n, iters);
prop_prob_nets = NaN(n, n, iters);
prop_prob_nets_weighted = NaN(n, n, iters);

for i = 1:iters
    fprintf('Iteration %d of %d.\n', i, iters);
    const_prob_nets(:, :, i) = constant_probability(n, p);
    const_prob_nets_weighted(:, :, i) = ...
        make_weighted_from_order(const_prob_nets(:, :, i), 1:n);
    prop_prob_nets(:, :, i) = constant_probability(n, p);
    prop_prob_nets_weighted(:, :, i) = ...
        make_weighted_from_order(prop_prob_nets(:, :, i), 1:n);
end

save_string = fullfile(base_path, ...
    'v7/Data/Simulations/Preprocessed/simulated_nets.mat');
save(save_string, 'const_prob_nets', 'const_prob_nets_weighted', ...
    'prop_prob_nets', 'prop_prob_nets_weighted');