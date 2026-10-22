% Predicted Centroid and Close Pair
% Mikko I. Malinen
% AI has been used to create the subalgorithms except the proposed.
% 2026
% v0.02

clear all;

ALGORITHM = '5'; % '1a' = centrality_filter_closest_pair 
                  % '1b' = geometric_median_pair
                  % '1c' = graph_based_pair
                  % '2' = weighted_scalarization
                  % '3' = constraint_based
                  % '4' = pareto_front_pairs
                  % '5' = Proposed

% The following line affects only the Algorithm '5' 
METHOD = 1;   % 1 = multidim by complex numbers (Proposed), 
              % 2 = Euclidean (Matlab), 
              % 3 = fastEuclidean (Matlab)  
              % 4 = fast 2-d Euclidean from Google Gemini/Stack Overflow site

%  data points

% X = [0 0; 0 2; 5 0; 5 1];
%X = load('../../datasets/s2.txt');

%X = randi(1000,10,2);    % range = 1..1000  n = 10, d = 2 

%X = randi(1000,100,2);    % range = 1..1000  n = 100, d = 2 

%X = randi(1000,100,3);    % range = 1..1000  n = 100, d = 3 

% 1000 points
%X = randi(1000,1000,2);    % range = 1..1000  n = 1000, d = 2 

%X = load('../../datasets/s2.txt');
%X = load('../../datasets/s4.txt');
%X = load('../../datasets/unbalance.txt');
%X = load('../../datasets/a3.txt');
%X = load('../../datasets/users_locations_Joensuu.txt');
%X = load('../../datasets/users_locations_Finland_until_2012.txt');

% Multi-dimensional datasets

%X = load('../../datasets/DIM032.txt'); % 32-dim, n= 1024
X = load('../../datasets/wine.txt'); % n=178, d=13
%X = load('../../datasets/bridge.txt'); % n=4096, d=16
%X = load('../../datasets/letter.txt'); % n=20000, d=16

n = size(X,1);
d = size(X,2);

if ALGORITHM == '5'

tic;

if METHOD == 1      % Transform into complex numbers
    if d/2 == floor(d/2)   % even
        ii = 1:n;
        d1 = 2:2:d;   % for poistettu
            z(ii,d1/2) = complex(X(ii,d1-1),X(ii,d1));
        %end
    else
        % Add a zero column to X
        X2 = zeros(n,d+1);
        X2(:,1:d) = X;
        ii = 1:n;
        d1 = 2:2:d+1  % for removed
            z(ii,d1/2) = complex(X2(ii,d1-1),X2(ii,d1));
        %end
    end
end % METHOD == 1


% 4 points
%
% z = [complex(0,0); complex(0,2); complex(5,0); complex(5,1)];  % first, (1,2) ok
%z = [complex(0,0); complex(0,10); complex(5,5); complex(6,1)]; % (1,2) ok
%z = [complex(0,0); complex(10,10); complex(5,5); complex(6,1)]; % (2,3) ok

%6 points
%z = [complex(0,0); complex(0,1); complex(5,0); complex(5,1); complex(6,3); complex(2,2)];

for repeats = 1:1      % repeats may be set higher for running time
                       % calculation, but every repeat gives same result

%L = X(:,:)*X(:,:)';

%size(L)

if METHOD == 1
    if d/2 == floor(d/2)   % even number of dimensions
        L = zeros(n,n);
        for d1 = 1:(d/2)
            ii = 1:n;
            L = L +  abs(z(ii,d1)- z(ii,d1).').^2;
        end
    else % odd number of dimensions
        L = zeros(n,n);
        for d1 = 1:((d+1)/2)
            ii = 1:n;
            L = L +  abs(z(ii,d1)- z(ii,d1).').^2;
        end
    end
    L = sqrt(L);  % This line can be removed without affecting the order of points
    %L = abs(z - z.');
elseif METHOD == 2
D = pdist(X);
L = squareform(D);
elseif METHOD == 3
D = pdist(X,'fasteuclidean',CacheSize=1e3);
L = squareform(D);
elseif METHOD == 4

% X is nx2
% Compute the squared norms of each row
rowNorms = sum(X.^2, 2); 

% Use implicit broadcasting to calculate the squared distance matrix
% (Resulting matrix will be nxn)
L_squared = rowNorms + rowNorms' - 2 * (X * X');

% Clean up negative floating-point zeroes, then square root
L = sqrt(max(L_squared, 0));

toc;

% Convert to a 1x12497500 row vector if you need to match pdist output format
%L = L(triu(true(size(L)), 1))';

end; 


L = rescale(L); % rescale to [0,1]

%-----------
% --- Setup Mock Data for a 10,000 x 10,000 Matrix ---
N = n; %10000;% Simulating a scaled distance matrix with some near-0 and near-1 values
%L = rand(N, N); 
% --- Step 1: Handle Overflows (Stabilize the Boundaries) ---% Define an epsilon floor to prevent log(0) -> Inf
eps_floor = 1e-15;  % 1e-15; 
L(L < eps_floor) = eps_floor; 
% --- Step 2: Handle Underflow/Precision Loss Near 1 ---% For values very close to 1, standard log loses precision.% We calculate the delta from 1 and use log1p for maximum accuracy.
Delta = L - 1;
% --- Step 3: Vectorized Information Calculation ---% log1p(x) calculates log(1 + x). Since Delta = L - 1, % log1p(L - 1) is mathematically identical to log(L), but precision-safe.
InfoMatrix = -log1p(Delta);
% Force the diagonal (self-distance) to 0 so a node doesn't measure info to itself
InfoMatrix(1:N+1:end) = 0; 
% --- Step 4: Precision-Safe Aggregation (Pairwise Summation) ---% DO NOT use a 'for' loop here. MATLAB's native sum() uses a pairwise % cascade algorithm under the hood which completely mitigates accumulator drift.
total_info_per_node = sum(InfoMatrix, 2); 
% --- Step 5: Identify the Minimum Information Node ---
[min_info_value, min_info_node_index] = min(total_info_per_node);
[max_info_value, max_info_node_index] = max(total_info_per_node);
[max_info_value, max_info_node_index] = max(total_info_per_node);
[B, info_node_indices_sorted] = sort(total_info_per_node,'descend');

toc;

fprintf('\nThe minimal node is: %d\n', min_info_node_index);
fprintf('Total Information: %f\n\n', min_info_value);
fprintf('The maximal node is: %d\n', max_info_node_index);
fprintf('Total Information: %f\n\n', max_info_value);

%--------------
figure;
scatter(X(:,1),X(:,2),'b.')
%scatter3(X(:,1),X(:,2),X(:,3),'b.')
hold on;
%sz = 100:-floor(100/n):floor(100/n);
%sz = B./B(1).*100;
scatter(X(min_info_node_index,1),X(min_info_node_index,2),100,'bO','filled')   % size sz shows centrality and closeness
%scatter3(X(indices,1),X(indices,2),X(indices,3),sz,'bO','filled')   % size sz shows centrality and closeness
%title('Dataset')
hold on;
scatter(X(info_node_indices_sorted(2),1),X(info_node_indices_sorted(2),2),100,[0.9290 0.6940 0.1250],'filled')% 100 is size
hold on;
scatter(X(max_info_node_index,1),X(max_info_node_index,2),100,'rO','filled')% 100 is size
%hold on;
%scatter(X(info_node_indices_sorted(2),1),X(info_node_indices_sorted(2),2),100,[0.9290 0.6940 0.1250],'filled')% 100 is size
%scatter3(X(indices(1),1),X(indices(1),2),X(indices(1),3),100,'rO','filled')       % 100 is size
%hold on;
% 100-floor(100/n) is size, [0.9290 0.6940 0.1250] is orange
%max_info_node_indexscatter(X(indices(2),1),X(indices(2),2),100-floor(100/n),[0.9290 0.6940 0.1250],'filled')     
%scatter3(X(indices(2),1),X(indices(2),2),X(indices(2),3),100-floor(100/n),[0.9290 0.6940 0.1250],'filled')     
title('Predicted Centroid (red), Close pair (red and orange) and Outlier (blue)');


end % repeats



elseif ALGORITHM == '1a'

    
    k = 10;
    [p,q] = centrality_filter_closest_pair(X,k)

elseif ALGORITHM == '1b';

    k = 10;
    [p,q] = geometric_median_pair(X,k)

elseif ALGORITHM == '1c'

    k = 10;
    knn_k = 4;
    [p,q] = graph_based_pair(X,k,knn_k)

elseif ALGORITHM == '2'

    
    lambda = 0.5;
    [p,q] = weighted_scalarization(X,lambda)
   

elseif ALGORITHM == '3'

    tau = (max(X(:,1))-min(X(:,1)))/10; % Possibly change 
    [p,q] = constraint_based(X,tau);

elseif ALGORITHM == '4'

    
    pareto_pairs = pareto_front_pairs(X)
    

end % ALGORITHM


