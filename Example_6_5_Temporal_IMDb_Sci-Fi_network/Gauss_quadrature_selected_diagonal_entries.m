% Description: computes upper Gauss--Radau quadrature bounds for a
% specified number of node-layer pairs top-ranked in diagonal
% estimation (the file diagonal_estimation_SC_Hadamard.m in the same 
% directory must be run first to have the diagonal estimation result 
% available) for the temporal IMDb Sci-Fi network. The accurate 
% approximations returned by Gauss--Radau quadrature are used to evaluate 
% the errors of diagonal estimation that are displayed in [1, Fig. 10]. 
% Note that, as the diagonal estimator yields upper bounds for non-negative
% matrices the comparison of both approximations allows the reliable
% identification of the top-ranked node-layer pairs.
%
% [1] Bergermann, K. & Stoll, M. (2021) Fast computation of matrix 
% function-based centrality measures for layer-coupled multiplex networks. 
% https://arxiv.org/abs/2104.14368v3
% 
% Kai Bergermann, 2021

%% load subroutines
addpath('../subroutines')
addpath('../subroutines/funm_kryl')

%% load results from diagonal estimation
load results/diagonal_estimation_SC.mat

% decive for broadcaster or receiver
broadcaster = 1;

% sort corresponding diagonal estimations
if broadcaster==1
    [D_sorted, D_sorted_ind] = sort(D_broadcaster, 'descend');
else
    [D_sorted, D_sorted_ind] = sort(D_receiver, 'descend');
end

% specify the number of top-ranked nodes in diagonal estimation for which
% to compute Gauss quadrature bounds
n_nodes=100;

gauss_quadrature_results = zeros(n_nodes,1);

%% load intra-layer adjacency matrices
load data/sci_fi_seasons_intra_layer_adjacency.mat

% extract network size
L=size(layer_id_table,1);
n=size(node_id_table,1);
nL=n*L;

%% assemble supra-adjacency matrix
deltaT = layer_id_table(2:end,2)-layer_id_table(1:end-1,2);
inter_layer_adjacency_matrix = spdiags([zeros(L,1),[exp(-deltaT);1]],0:1,L,L);

omega=10^1;

A = A_intra + omega * kron(inter_layer_adjacency_matrix,speye(n));

% spy plot of supra-adjacency matrix
figure(1)
spy(A)

%% estimating the extremal eigevalues via an initial lanczos procedure
beta=1e-02;

% set funm_kryl parameters
param.function = @(A) expm(beta*A);
param.restart_length = 5;
param.max_restarts = 50;
param.hermitian = 0;          % set 0 if A is not Hermitian
param.V_full = 0;             % set 1 if you need Krylov basis
param.H_full = 1;             % if using rational functions you can set this 0
param.exact = [];          % if not known set to []
param.bound = 1;              % returns upper and lower bounds (after some cycles)
param.stopping_accuracy = 1e-12;  % stopping accuracy
param.inner_product = @inner_product;
param.thick = [];             % thick-restart function  
param.min_decay = 0.95;       % we desire linear error reduction of rate < .95 
param.waitbar = 1;            % show waitbar 
param.reorth_number = 0;      % #reorthogonalizations
param = param_init(param);    % check and correct param structure

[~,out] = funm_kryl(A,ones(nL,1),param);
[~,lambda]=eig(out.H_full);

lambda_max=max(max(lambda));
lambda_min=min(min(lambda));

%% assemble symmetric bipartite network representation
A_bipartite=[sparse(nL,nL), A; A', sparse(nL,nL)];

% set parameter beta
beta = 5/lambda_max;

% set number of Lanczos iterations
maxit=10;

%% loop over selected node-layer pairs
for k=1:n_nodes
    tic
    
    % set unit vectors as starting vectors
    ei=zeros(2*nL,1); ei(D_sorted_ind(k))=1;
    % Lanczos tridiagonalization
    T_jp1=lanczos_tridiag(A_bipartite,ei,maxit);
    
    % Compute lower Gauss--Radau bound on e_i^T*f(A)*e_i
    gauss_quadrature_results(k) = gauss_radau_subgraph(T_jp1,beta,lambda_min);
    toc
end

% output array with Gauss quadrature results for the top 'n_nodes' joint
% centralities from diagonal estimation
results = [gauss_quadrature_results, D_sorted(1:n_nodes), D_sorted_ind(1:n_nodes)]
