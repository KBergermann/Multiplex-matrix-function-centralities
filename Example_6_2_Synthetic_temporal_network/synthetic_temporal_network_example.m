% Description: This script creates a small synthetic temporal network with 
% n=200 nodes and L=4 time layers according to the description in 
% [1, Sec. 5.2]. The random edges of one 'agenda-setting' node in the
% first time layer are replaced by paths that distribute information
% through the network more effectively than random edges would. 
% This script reproduces [2, Fig. 6] in which we compare dynamic 
% communicabilities with marginal node Katz centralities for forward and
% backward evolving time, plotted against the aggregate node degree.
%
% [1] Fenu, C. & Higham, D.J. (2017) Block matrix formulations for evolving networks,
% SIAM Journal on Matrix Analysis and Applications, vol. 38, no. 2, pp. 343-360.
% https://doi.org/10.1137/16M1076988
% [2] Bergermann, K. & Stoll, M. (2021) Fast computation of matrix 
% function-based centrality measures for layer-coupled multiplex networks. 
% https://arxiv.org/abs/2104.14368v3
% 
% Kai Bergermann, 2021

%% create network
% create the random network as described in [1, Sec. 5.2]
n=200;
L=4;

% create directed random graph with probability 4/n for the presence of an edge
rng(28) 
R = rand(n,n,L);
A = double(R<4/n);
% remove all edges that emanate from node 1
A(1,:,:) = 0;

% create 16 paths through time starting at node 1
for i=1:16
    n2 = randi([1,n]);
    A(1,n2,1) = 1;
    n3 = randi([1,n]);
    A(n2,n3,2) = 1;
    n4 = randi([1,n]);
    A(n3,n4,3) = 1;
    n5 = randi([1,n]);
    A(n4,n5,4) = 1;
end

% delete repeated edges within a time level and delete self-loops
A(A>1) = 1;
A(:,:,1) = A(:,:,1) - diag(diag(A(:,:,1)));
A(:,:,2) = A(:,:,2) - diag(diag(A(:,:,2)));
A(:,:,3) = A(:,:,3) - diag(diag(A(:,:,3)));
A(:,:,4) = A(:,:,4) - diag(diag(A(:,:,4)));

%% compute dynamic centralities

% aggregate out-degree
aggregate_degree = sum(sum(A,2),3);

% determine maximum spectral radius over the time levels
lambda = zeros(L,1);
for l=1:L
    lambda(l) = eigs(A(:,:,l),1,'lm');
end
p_star = max(lambda);

% compute dynamic broadcast communicability
b = 1;
deltaT = 1;
alpha = 0.9/p_star;
Sj = zeros(n,n);
for l=1:L
    Sj = (eye(n,n) + exp(-b*deltaT)*Sj)*((eye(n,n) - alpha*A(:,:,l))\eye(n,n)) - eye(n,n);
end
dynamic_broadcast_communicability = sum(Sj/norm(Sj,2),2);

% plot aggregate_degree vs. dynamic_broadcast_communicability 'Original'
% (forward time)
figure(1)
subplot(221)
scatter(aggregate_degree(2:end),dynamic_broadcast_communicability(2:end),'b.')
hold on
scatter(aggregate_degree(1),dynamic_broadcast_communicability(1),'r*')
title('Original')
xlabel('Aggregate degree')
ylabel('Dynamic communicability')
ylim([0 inf])
hold off

% compute reverse dynamic broadcast communicability
b = 1;
deltaT = 1;
alpha = 0.9/p_star;
Sj = zeros(n,n);
for l=L:-1:1
    Sj = (eye(n,n) + exp(-b*deltaT)*Sj)*((eye(n,n) - alpha*A(:,:,l))\eye(n,n)) - eye(n,n);
end
reverse_dynamic_broadcast_communicability = sum(Sj/norm(Sj,2),2);

% plot aggregate_degree vs. reverse_dynamic_broadcast_communicability
% 'Time-reversed' (backward time)
subplot(222)
scatter(aggregate_degree(2:end),reverse_dynamic_broadcast_communicability(2:end),'b.')
hold on
scatter(aggregate_degree(1),reverse_dynamic_broadcast_communicability(1),'r*')
title('Time-reversed')
xlabel('Aggregate degree')
ylabel('Dynamic communicability')
ylim([0 inf])
hold off


%% compute multiplex Katz centrality

% set interlayer adjacency
interlayer_type=2;

switch interlayer_type
    case 1 % temporal (bi-directional) -- reproduces lower two plots of  [1, Fig. 2]
        interlayer_adjacency_matrix = spdiags([ones(L,1),zeros(L,1),ones(L,1)],-1:1,L,L); 
    case 2 % temporal (uni-directional [forward])
        interlayer_adjacency_matrix = spdiags([zeros(L,1),ones(L,1)],0:1,L,L);
end

% build supra-adjacency matrix with time layer ordering 1, 2, 3, 4
A_supra = blkdiag(A(:,:,1), A(:,:,2), A(:,:,3), A(:,:,4)) + kron(interlayer_adjacency_matrix, eye(n,n));
nL = n*L;

% compute largest eigevalue and set parameter alpha
lambda_supra = eigs(A_supra,1,'lm');
alpha_supra = 0.9/lambda_supra;

% compute marginal node Katz centrality 'Original'
supracentrality_original = sum(reshape(sum((eye(nL,nL) - alpha_supra*A_supra)\eye(nL,nL),2),n,L),2);

% plot aggregate_degree vs. supracentrality_original 'Original'
% (forward time)
subplot(223)
scatter(aggregate_degree(2:end),supracentrality_original(2:end),'b.')
hold on
scatter(aggregate_degree(1),supracentrality_original(1),'r*')
title('Original')
xlabel('Aggregate degree')
ylabel('Katz MNC')
ylim([0 inf])
hold off

% build supra-adjacency matrix with time layer ordering 4, 3, 2, 1
A_supra_reversed = blkdiag(A(:,:,4), A(:,:,3), A(:,:,2), A(:,:,1)) + kron(interlayer_adjacency_matrix, eye(n,n));

% compute largest eigevalue and set parameter alpha
lambda_supra_reversed = eigs(A_supra_reversed,1,'lm');
alpha_supra_reversed = 0.9/lambda_supra_reversed;

% compute marginal node Katz centrality 'Time-reversed'
supracentrality_reversed = sum(reshape(sum((eye(nL,nL) - alpha_supra_reversed*A_supra_reversed)\eye(nL,nL),2),n,L),2);

% plot aggregate_degree vs. supracentrality_reversed 'Time-reversed'
% (backward time)
subplot(224)
scatter(aggregate_degree(2:end),supracentrality_reversed(2:end),'b.')
hold on
scatter(aggregate_degree(1),supracentrality_reversed(1),'r*')
title('Time-reversed')
xlabel('Aggregate degree')
ylabel('Katz MNC')
ylim([0 inf])
hold off
