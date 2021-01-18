function is_connected = is_fully_connected(G)
% is_fully_connected Determines if the network given by the adjacency
% matrix G is fully connected.

g = digraph(G);
bins = conncomp(g, 'Type', 'Weak');
is_connected = all(bins == 1);

end

