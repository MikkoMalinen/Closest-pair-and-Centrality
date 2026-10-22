function paretoPairs = pareto_front_pairs(X)
% PARETO FRONT ENUMERATION
% References: Multiobjective Pareto theory [hashem2025rose]

tic;

n = size(X,1);

% Centrality: average distance
C = zeros(n,1);
for i = 1:n
    C(i) = mean(vecnorm(X - X(i,:), 2, 2));
end

pairs = [];
for i = 1:n
    for j = i+1:n
        d = norm(X(i,:) - X(j,:));
        c = C(i) + C(j);
        pairs = [pairs; i, j, d, c];
    end
end

% Pareto filtering
pareto = [];
for a = 1:size(pairs,1)
    dominated = false;
    for b = 1:size(pairs,1)
        if b ~= a
            if pairs(b,3) <= pairs(a,3) && pairs(b,4) >= pairs(a,4) ...
               && (pairs(b,3) < pairs(a,3) || pairs(b,4) > pairs(a,4))
                dominated = true;
                break;
            end
        end
    end
    if ~dominated
        pareto = [pareto; pairs(a,:)];
    end
end

toc;

figure;
scatter(X(:,1),X(:,2),10,'b.')
hold on;
scatter(X(pareto(:,1),1),X(pareto(:,1),2),20,'rO','filled')
hold on;
scatter(X(pareto(:,2),1),X(pareto(:,2),2),20,[0.9290 0.6940 0.1250],'filled')


paretoPairs = pareto;
end
