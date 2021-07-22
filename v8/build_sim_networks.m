clc; close all; clear;

%% Specify paths and miscellaneous settings

base_path = '/Volumes/My Passport/Curiosity/';
addpath(genpath(fullfile(base_path, 'Helper')))

%% Build networks

n = 70; % number of nodes

p = 0.4; % parameter for constant probability model
m_0 = 4; % parameter for preferential attachment model
m = 4; % parameter for preferential attachment model

iters = 100;

const_prob_nets = NaN(n, n, iters);
const_prob_nets_weighted = NaN(n, n, iters);
prop_prob_nets = NaN(n, n, iters);
prop_prob_nets_weighted = NaN(n, n, iters);
PA_nets = NaN(n, n, iters);
PA_nets_weighted = NaN(n, n, iters);

for i = 1:iters
    fprintf('Iteration %d of %d.\n', i, iters);
    const_prob_nets(:, :, i) = constant_probability(n, p);
    const_prob_nets_weighted(:, :, i) = ...
        make_weighted_from_order(const_prob_nets(:, :, i), 1:n);
    prop_prob_nets(:, :, i) = proportional_probability(n);
    prop_prob_nets_weighted(:, :, i) = ...
        make_weighted_from_order(prop_prob_nets(:, :, i), 1:n);
    PA_nets(:, :, i) = preferential_attachment(n, m_0, m);
    PA_nets_weighted(:, :, i) = ...
        make_weighted_from_order(PA_nets(:, :, i), 1:n);
end

save_string = fullfile(base_path, ...
    'v8/Data/Simulations/simulated_nets.mat');
save(save_string, 'const_prob_nets', 'const_prob_nets_weighted', ...
    'prop_prob_nets', 'prop_prob_nets_weighted', ...
    'PA_nets', 'PA_nets_weighted');