function [p, q] = graph_based_pair(X, k, knn_k)
% GRAPH-BASED SEARCH USING kNN GRAPH
% References: Graph-based clustering [li2012graphclustering],
%             HNSW proximity graphs [malkov2018hnsw]

tic;

n = size(X,1);

% Build kNN graph
[idx, ~] = knnsearch(X, X, 'K', knn_k+1);
idx = idx(:,2:end); % remove self

% Graph centrality: degree (local density)
degree = sum(idx > 0, 2);

% Select top k central nodes
[~, topIdx] = maxk(degree, k);
S = X(topIdx,:);

% Find closest adjacent pair
bestDist = inf;
p = []; q = [];

for i = 1:k
    xi = S(i,:);
    for j = 1:k
        if i ~= j
            d = norm(xi - S(j,:));
            if d < bestDist
                bestDist = d;
                p = xi;
                q = S(j,:);
            end
        end
    end
end

toc;

end

