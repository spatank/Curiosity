function [G, node_order] = proportional_probability(n)
% Probability of edge formation proportional to node number added
%   n: network size (number of nodes)
G = zeros(n);
for i = 1:n
    for j = 1:i-1
        r = rand;
        if r < i/n
            G(i, j) = 1;
            G(j, i) = 1;
        end
    end
end
node_order = 1:n;
end

