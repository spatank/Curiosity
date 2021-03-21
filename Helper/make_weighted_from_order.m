function weighted_G = make_weighted_from_order(G, node_order)
% Adapted from code by ASB
%
% Original at https://github.com/BassettLab/Reorderability_scripts
% This function returns a matrix with ranks assigned to network edges.
% 
% Edges between nodes that are 'born' earlier have higher ranks 
% (i.e. smaller values). Non-existent edges take the worst rank possible 
% (i.e. Inf).
%
% Arguments:
% G: adjacency matrix
% node_order: ranks assigned to the nodes (usually time ordering)
%
% Outputs:
% weighted_G: rank-ordered matrix

reordered_G = G(node_order, node_order); % often unnecessary
n = length(node_order); % number of nodes
val_mat = ones(n, n);
for col = 1:n
    val_mat(1:col, col) = val_mat(1:col, col) * col;
    val_mat(col, 1:col) = val_mat(col, 1:col) * col;
end
weighted_G = reordered_G .* val_mat;
weighted_G(logical(eye(n))) = node_order; 
% replace 0 weighted edges with the largest edge weight possible
% this is equivalent to assigning these edges the worst rank possible
weighted_G(weighted_G == 0) = Inf;


end

