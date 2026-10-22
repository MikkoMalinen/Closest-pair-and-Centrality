% Tentative Centroid and Closest Pair by real (float) numbers
% Mikko I. Malinen
% 2026
% v0.01

clear all;



%  data points

% X = [0 0; 0 2; 5 0; 5 1];
%X = load('../../datasets/s2.txt');

%X = randi(1000,10,2);    % range = 1..1000  n = 10, d = 2 
X = randi(1000,100,2);    % range = 1..1000  n = 100, d = 2 
% 1000 points is too much:
%X = randi(1000,1000,2);    % range = 1..1000  n = 1000, d = 2 

% Let's transform X to complex numbers

n = size(X,1);

%z = [complex(X(1,1),X(1,2)); complex(X(2,1),X(2,2)); complex(X(3,1),X(3,2)); complex(X(4,1),X(4,2))];

for ii = 1:n
    z(ii) = complex(X(ii,1),X(ii,2));
end


% 4 points
%
% z = [complex(0,0); complex(0,2); complex(5,0); complex(5,1)];  % first, (1,2) ok
%z = [complex(0,0); complex(0,10); complex(5,5); complex(6,1)]; % (1,2) ok
%z = [complex(0,0); complex(10,10); complex(5,5); complex(6,1)]; % (2,3) ok

%6 points
%z = [complex(0,0); complex(0,1); complex(5,0); complex(5,1); complex(6,3); complex(2,2)];

for repeats = 1:1      % repeats may be set higher for running time
                       % calculation, but every repeat gives same result

tic;


number_of_clusters = 2;


cluster_size = floor(n/number_of_clusters);
for k = 1:n  % k is here index of point inside contour
product = ones(n,1);
for point_number = 1:n
    if k ~= point_number
        product(k) = product(k) * sqrt((X(k,:)-X(point_number,:))*(X(k,:)-X(point_number,:))');
    end
end
Res(k) = 1/product(k);
end % for k

toc;

end % repeats

%for ii = 1:n
%    Abs_of_Res(ii) = abs(Res(ii));
%end
Abs_of_Res = Res;

[B,indices] = sort(Abs_of_Res,'descend')

figure;
scatter(X(:,1),X(:,2),'b.')
hold on;
%sz = 100:-floor(100/n):floor(100/n);
sz = B./B(1).*100;
scatter(X(indices,1),X(indices,2),sz,'bO','filled')   % size sz shows centrality and closeness
%title('Dataset')
hold on;
scatter(X(indices(1),1),X(indices(1),2),100,'rO','filled')       % 100 is size
hold on;
% 100-floor(100/n) is size, [0.9290 0.6940 0.1250] is orange
scatter(X(indices(2),1),X(indices(2),2),100-floor(100/n),[0.9290 0.6940 0.1250],'filled')     
title('Tentative Centroid (red) and Closest pair (red and orange)');
