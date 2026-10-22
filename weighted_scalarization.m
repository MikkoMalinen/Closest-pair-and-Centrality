function [p, q] = weighted_scalarization(X, lambda)
% WEIGHTED SCALARIZATION
% J = lambda * ||p - q|| + (1-lambda) * (-C(p) - C(q))
% References: Multiobjective scalarization [hashem2025rose]

tic;

n = size(X,1);

% Centrality: average distance
C = zeros(n,1);
for i = 1:n
    C(i) = mean(vecnorm(X - X(i,:), 2, 2));
end

bestVal = inf;
p = []; q = [];

for i = 1:n
    for j = i+1:n
        d = norm(X(i,:) - X(j,:));
        J = lambda * d + (1-lambda) * (C(i) + C(j));
        if J < bestVal
            bestVal = J;
            p = X(i,:);
            q = X(j,:);
        end
    end
end

toc;

figure;
scatter(X(:,1),X(:,2),10,'b.')
hold on;
scatter(p(1),p(2),20,'rO','filled')
hold on;
scatter(q(1),q(2),20,[0.9290 0.6940 0.1250],'filled')

end
