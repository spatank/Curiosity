clc; close all; clear;

load('/Volumes/My Passport/Curiosity/v4/Data/KNOT/Raw/subj_101.mat');

adj = double(adj);
n = size(adj, 1);
[components, component_sizes] = conncomp(digraph(adj), 'Type', 'Weak');

length(unique(components))

iters = 100;
nulls = zeros(n, n, iters);

i = 1;
counter = 1;
while i < iters
    candidate_net = randmio_und(adj, 10);
    [components, component_sizes] = conncomp(digraph(candidate_net), 'Type', 'Weak');
    if length(unique(components)) == 1
        nulls(:, :, i) = candidate_net;
        i = i + 1;
    end
    counter = counter + 1;
end
