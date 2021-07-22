function [G, node_order] = preferential_attachment(n, m_0, m)
% Barabasi-Albert preferential attachment model
%   n: network size (number of nodes)
%   m_0: size of connected random network to start with
%   m: number of nodes to preferentially connect to at each step
G = zeros(n);
% construct random network
is_connected = 0;
while is_connected ~= 1
    rand_net = rand(m_0) - eye(m_0);
    rand_net = (rand_net + rand_net')/2;
    rand_net(rand_net < 0.5) = 0;
    rand_net(rand_net > 0) = 1;
    % ensure network is connected
    [components, ~] = conncomp(digraph(rand_net));
    is_connected = max(components);
end
G(1:m_0, 1:m_0) = rand_net;
% add nodes 
for n = (m_0 + 1):n
    P = sum(G)./sum(G(:));
    cP = cumsum(P);
    for m_i = 1:m
        % choose node i
        r = rand;
        node_i = find(r < cP, 1, 'first');
        while ismember(node_i, find(G(n, :)))
            r = rand;
            cP = cumsum(P);
            node_i = find(r < cP, 1, 'first');
        end
        G(n, node_i) = 1;
        G(node_i ,n) = 1;
    end
end
node_order = 1:n;
end