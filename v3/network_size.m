function [num_nodes, num_edges] = network_size(adj)
adj = double(adj);
adj(adj == 0) = NaN;
edg = Adj2Edg(adj);
num_nodes = size(adj, 1);
num_edges = length(edg)/2;
end

