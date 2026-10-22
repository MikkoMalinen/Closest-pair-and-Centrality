function [p, q] = centrality_filter_closest_pair(X, k)
% CENTRALITY FILTER + CLOSEST PAIR
% Centrality measure: average distance to all points
% References: Tukey depth [tukey1975depth], medoid centrality [rousseeuw1987medoid]

tic;

n = size(X,1);

% Compute centrality scores
C = zeros(n,1);
for i = 1:n
    C(i) = mean(vecnorm(X - X(i,:), 2, 2));
end

% Select k most central points (smallest average distance)
[~, idx] = mink(C, k);
S = X(idx,:);

% Compute closest pair inside S
[p, q] = closest_pair(S);

toc;

figure;
scatter(X(:,1),X(:,2),10,'b.')
hold on;
scatter(p(1),p(2),20,'rO','filled')
hold on;
scatter(q(1),q(2),20,[0.9290 0.6940 0.1250],'filled')

end


function [p, q] = closest_pair(S)
% Standard O(k^2) closest pair (sufficient for moderate k)
k = size(S,1);
bestDist = inf;
p = []; q = [];

for i = 1:k-1
    for j = i+1:k
        d = norm(S(i,:) - S(j,:));
        if d < bestDist
            bestDist = d;
            p = S(i,:);
            q = S(j,:);
        end
    end
end
end
