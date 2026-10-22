function [p, q] = constraint_based(X, tau)
% CONSTRAINT-BASED OPTIMIZATION
% Minimize ||p - q|| subject to C(p) >= tau/2,   C(q) >= tau/2
% References: Depth-based centrality [tukey1975depth]

tic;

n = size(X,1);

% Centrality: average distance
C = zeros(n,1);
for i = 1:n
    %C(i) = mean(vecnorm(X - X(i,:), 2, 2));
    C(i) = vecnorm(mean(X) - X(i,:),2,2);
end

% Filter points
idx = find(C <= tau/2);   % smaller C = more central
S = X(idx,:);

% Closest pair in filtered set
[p, q] = closest_pair(idx,S);

toc;

figure;
scatter(X(:,1),X(:,2),10,'b.')
hold on;
p = p
q = q
scatter(X(p,1),X(p,2),20,'rO','filled')
hold on;
scatter(X(q,1),X(q,2),20,[0.9290 0.6940 0.1250],'filled')

end

%function [p, q] = closest_pair(idx,S)
%distmin = Inf; p=0; q= 0;
%for i=1:length(idx)-1
%    for j=i+1:length(idx)
%        dist = (S(i,1)-S(j,1))^2+(S(i,2)-S(j,2))^2;
%        if dist <= distmin
%            p = idx(i);
%            q = idx(j);
%            distmin = dist;
%        end
%    end
%end
%end

function [p, q] = closest_pair(idx,S)
% Standard O(k^2) closest pair (sufficient for moderate k)
k = size(S,1);
bestDist = Inf;
p = []; q = [];

for i = 1:k-1
    for j = i+1:k
        d = norm(S(i,:) - S(j,:));
        if d < bestDist
            bestDist = d;
            p = idx(i);  % S(i,:);
            q = idx(j);  % S(j,:);
        end
    end
end
end