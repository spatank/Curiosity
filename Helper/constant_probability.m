function [G, node_order] = constant_probability(n, p)
% Edges form between node pairs with constant probability p
%   n: network size (number of nodes)
%   p: probability of edge formation between nodes
G = zeros(n);
for i = 1:n
    for j = 1:i-1
        r = rand;
        if r < p
            G(i, j) = 1;
            G(j, i) = 1;
        end
    end
end
node_order = 1:n;
end

