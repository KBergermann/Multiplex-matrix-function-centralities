% Description: computes lower and upper bounds on the Estrada index of the
% European airlines multiplex network [1,2] by applying Gauss,
% Gauss--Radau, and Gauss--Lobatto quadrature rules to all diagonal entries
% of f(A)=exp(beta*A) where A is the supra-adjacency matrix of the
% multiplex network. The bounds for 1 up to 5 Lanczos iterations are
% computed and displayed in [3, Tab. 2].
% 
% [1] Cardillo, A., Gómez-Gardenes, J., Zanin, M., Romance, M., Papo, D.,
% Del Pozo, F. & Boccaletti, S. (2013) Emergence of network features from 
% multiplexity. Sci. Rep., 3(1), 1-6.
% [2] Taylor, D. (2021) Code Release: Supracentrality. Available at 
% https://github.com/taylordr/Supracentrality.
% [3] Bergermann, K. & Stoll, M. (2021) Fast computation of matrix 
% function-based centrality measures for layer-coupled multiplex networks. 
% https://arxiv.org/abs/2104.14368v3
% 
% Kai Bergermann, 2021

%% load subroutines
addpath('../subroutines')
addpath('../subroutines/funm_kryl')

%% load intra-layer adjacency matrices as well as node and layer IDs
load data/multiplex_airlines_GC.mat
load data/multiplex_airlines_layer_ids.mat
load data/multiplex_airlines_airport_names_coords.mat

% extract network size
n=net.N;
L=net.T;
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
A_intra = spalloc(nL,nL,nL);
for t=1:L
  ids = (t-1)*n + (1:n);
  A_intra(ids,ids) = net.A{t}; %intra-layer block matrix
end

omegas = 10^(0);

A = A_intra + omegas(1) * kron(interlayer_adjacency_matrix,speye(n));

% spy plot of supra-adjacency matrix
figure(1)
spy(A)

%% visualization of layer l as Matlab graph
l=1;
G = graph(net.A{l}, 'OmitSelfLoops');
G.Nodes.Name=airport_names;
figure(2)
plot(G, 'LineWidth', 1, 'EdgeColor', [0, 0, 0],'XData',coords(:,1),'YData',coords(:,2));
title(layer_ids{l})

%% estimating the extremal eigevalues via an initial lanczos procedure
beta=1e-02;

% set funm_kryl parameters
param.function = @(A) expm(beta*A);        % other choices: 'expBA', 'expCF', ...
param.restart_length = 10;
param.max_restarts = 50;
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

[~,out] = funm_kryl(A,ones(nL,1),param);
lambda=eig(out.H_full);

lambda_min=min(lambda);
lambda_max=max(lambda);

%% compute subgraph centrality of all node-layer pairs
% set parameter beta
beta=5/lambda_max;

% specify range of iterations
maxit=1:5;
lower_gauss_subgraph=zeros(nL,length(maxit));
lower_gauss_radau_subgraph=zeros(nL,length(maxit));
upper_gauss_radau_subgraph=zeros(nL,length(maxit));
upper_gauss_lobatto_subgraph=zeros(nL,length(maxit));

% loop over range of iterations to approximate SC with Krylov subspace
% methods
for j=1:length(maxit)
    bar=waitbar(0,'Looping over nodes');

    tic
    for i=1:nL
        % set unit vectors as starting vectors
        ei=zeros(nL,1); ei(i)=1;
        % Lanczos tridiagonalization
        T_jp1=lanczos_tridiag(A,ei,maxit(j));
        
        % Compute lower and upper bounds on e_i^T*exp(-beta*(-A))*e_i
        lower_gauss_subgraph(i,j)=gauss_subgraph(T_jp1,beta);
        lower_gauss_radau_subgraph(i,j)=gauss_radau_subgraph(T_jp1,beta,lambda_min);
        upper_gauss_radau_subgraph(i,j)=gauss_radau_subgraph(T_jp1,beta,lambda_max);
        upper_gauss_lobatto_subgraph(i,j)=gauss_lobatto_subgraph(T_jp1,beta,lambda_min,lambda_max);
        waitbar(i/nL,bar);
    end
    close(bar)
    toc
end

%% print bounds on Estrada index as a function of the Lanczos iterations
fprintf('Bounds on the Estrada index:\n')
fprintf('| # iterations |      G (lower)      |     GR (lower)      |     GR (upper)      |     GL (upper)      |\n')
fprintf('|------------------------------------------------------------------------------------------------------|\n')

fprintf('|      %d       | %.0f               | %.0f               | %.0f               | %.0f              |\n',1,sum(lower_gauss_subgraph(:,1)),sum(lower_gauss_radau_subgraph(:,1)),sum(upper_gauss_radau_subgraph(:,1)),sum(upper_gauss_lobatto_subgraph(:,1)))
fprintf('|      %d       | %.5f         | %.5f         | %.5f         | %.5f         |\n',2,sum(lower_gauss_subgraph(:,2)),sum(lower_gauss_radau_subgraph(:,2)),sum(upper_gauss_radau_subgraph(:,2)),sum(upper_gauss_lobatto_subgraph(:,2)))
fprintf('|      %d       | %.9f     | %.9f     | %.9f     | %.9f     |\n',3,sum(lower_gauss_subgraph(:,3)),sum(lower_gauss_radau_subgraph(:,3)),sum(upper_gauss_radau_subgraph(:,3)),sum(upper_gauss_lobatto_subgraph(:,3)))
fprintf('|      %d       | %.12f  | %.12f  | %.12f  | %.12f  |\n',4,sum(lower_gauss_subgraph(:,4)),sum(lower_gauss_radau_subgraph(:,4)),sum(upper_gauss_radau_subgraph(:,4)),sum(upper_gauss_lobatto_subgraph(:,4)))
fprintf('|      %d       | %.12f  | %.12f  | %.12f  | %.12f  |\n',5,sum(lower_gauss_subgraph(:,5)),sum(lower_gauss_radau_subgraph(:,5)),sum(upper_gauss_radau_subgraph(:,5)),sum(upper_gauss_lobatto_subgraph(:,5)))
