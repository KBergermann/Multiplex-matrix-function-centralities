% Description: compares the approximation of total communicability obtained
% with Krylov subspace methods to the 'exact' quantity computed with expm
% for the weighted Scotland Yard network and plots the approximation error 
% used in [1, Fig. 7].
% 
% [1] Bergermann, K. & Stoll, M. (2021) Fast computation of matrix 
% function-based centrality measures for layer-coupled multiplex networks. 
% https://arxiv.org/abs/2104.14368v3
% 
% Kai Bergermann, 2021

%% load subroutines
addpath('../subroutines')
addpath('../subroutines/funm_kryl')

%% load intra-layer adjacency matrices
weighted_edges=1;
switch weighted_edges
    case 0
        load data/scotlandyard_adjacency.mat
        A_layer_adjacencies = A;
    case 1
        load data/scotlandyard_adjacency_weighted.mat
        A_layer_adjacencies=A_weighted;
end

% extract network size
L=size(A_layer_adjacencies,2);
n=size(A_layer_adjacencies{1},1);
nL=n*L;

%% choose inter-layer coupling
interlayer_type=1;
switch interlayer_type
    case 1 % all-to-all (self-edges excluded)
        interlayer_adjacency_matrix = sparse(ones(L)-eye(L));
    case 2 % all-to-all including self-edges
        interlayer_adjacency_matrix = ones(L);
    case 3 % temporal (uni-directional [forward])
        interlayer_adjacency_matrix = spdiags([zeros(L,1),ones(L,1)],0:1,L,L);
end

%% assemble supra-adjacency matrix
A_hat = spalloc(nL,nL,nL);
for t=1:L
  ids = (t-1)*n+(1:n);
  A_hat(ids,ids) = A_layer_adjacencies{t}; %intra-layer block matrix
end

omegas = 10^0;

A = A_hat + omegas(1)*kron(interlayer_adjacency_matrix,speye(n));

% spy plot of supra-adjacency matrix
figure(1)
spy(A)

% compute largest eigenvalue
lambda_max=eigs(A,1);

%% compute total communcability
% set parameter beta
beta=0.5/lambda_max;

% 'exact' result computed with expm
TC=expm(beta*A)*ones(nL,1);

% specify range of iterations
maxit=2:10;
TC_funm_kryl=zeros(nL,length(maxit));

% loop over range of iterations to approximate TC with Krylov subspace
% methods
for j=1:length(maxit)
    tic
    % set funm_kryl parameters
    param.function = @(A) expm(beta*A);        % other choices: 'expBA', 'expCF', ...
    param.restart_length = maxit(j);
    param.max_restarts = 1;
    param.hermitian = 1;          % set 0 if A is not Hermitian
    param.V_full = 0;             % set 1 if you need Krylov basis
    param.H_full = 1;             % if using rational functions you can set this 0
    param.exact = [];          % if not known set to []
    param.bound = 1;              % returns upper and lower bounds (after some cycles)
    param.stopping_accuracy = 1e-16;  % stopping accuracy
    param.inner_product = @inner_product;
    param.thick = [];             % thick-restart function  
    param.min_decay = 0.95;       % we desire linear error reduction of rate < .95 
    param.waitbar = 0;            % show waitbar 
    param.reorth_number = 0;      % #reorthogonalizations
    param = param_init(param);    % check and correct param structure

    % approximate TC according to [1, Sec. 5.1]
    [TC_funm_kryl(:,j),out] = funm_kryl(A,ones(nL,1),param);
    toc
end
    
%% error plot
figure(2)
semilogy(maxit,max(abs(TC-TC_funm_kryl)))

