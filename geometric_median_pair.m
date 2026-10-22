function [p, q] = geometric_median_pair(X, k)
% GEOMETRIC MEDIAN PAIR
% References: Robust centrality [rousseeuw1987medoid]

tic;

% Compute geometric median using Weiszfeld's algorithm
m = geometric_median(X);

% Sort points by distance to median
d = vecnorm(X - m, 2, 2);
[~, idx] = mink(d, k);
S = X(idx,:);

% Closest pair among S
[p, q] = closest_pair(idx,S);

toc;

figure;
scatter(X(:,1),X(:,2),10,'b.')
hold on;
scatter(X(p,1),X(p,2),20,'rO','filled')
hold on;
scatter(X(q,1),X(q,2),20,[0.9290 0.6940 0.1250],'filled')

end


function m = geometric_median(X)
% Weiszfeld's algorithm
m = mean(X,1);
for iter = 1:100
    d = vecnorm(X - m, 2, 2);
    w = 1 ./ max(d, 1e-12);
    m = sum(X .* w, 1) / sum(w);
end
end

function [p, q] = closest_pair(idx,S)
distmin = Inf; p=0; q= 0;
for i=1:length(idx)-1
    for j=i+1:length(idx)
        dist = (S(i,1)-S(j,1))^2+(S(i,2)-S(j,2))^2;
        if dist <= distmin
            p = idx(i);
            q = idx(j);
            distmin = dist;
        end
    end
end
end