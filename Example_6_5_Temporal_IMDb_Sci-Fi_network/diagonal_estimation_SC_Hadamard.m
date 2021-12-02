% Description: computes deterministic estimates of the diagonal of 
% f(A)=exp(beta*A) using Hadamard vetors in order to obtain approximations 
% to joint subgraph centralities of the temporal IMDb Sci-Fi network. 
% The results for different values of s (number of Hadamard vectors) are 
% compared with approximations computed to high precision with Gauss
% quadrature rules (see Gauss_quadrature_selected_diagonal_entries.m in the
% same directory) and displayed in [1, Fig. 10]. Furthermore, the estimates
% are used to generate the marginal node rankings in [1, Tab. 4].
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
load data/sci_fi_seasons_intra_layer_adjacency.mat

% extract network size
L=size(layer_id_table,1);
n=size(node_id_table,1);
nL=n*L;

%% assemble supra-adjacency matrix
% temporal inter-layer coupling with weights exp(-deltaT)
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

% update parameter beta in funm_kryl
param.function = @(A) expm(beta*A);
param.waitbar = 0;

%% subgraph centrality diagonal estimator
% using deterministic Hadamard vectors

tic
% build two small Hadamard matrices from which to compute larger ones via
% Kronecker products
H1 = [1, 1; 1, -1];
H2 = [1, 1, 1, 1; 1, -1, 1, -1; 1, 1, -1, -1; 1, -1, -1, 1];

% set number of Hadamard vectors
s=32;

% build Hadamard matrix of desired size by means of Kronecker products
% (H can be an arbitrary Hadamard matrix, we only need its orthogonality)
if s==2
    H = H1;
elseif s==4
    H = H2;
elseif s==8
    H = kron(H1,H2);
elseif s==16
    H = kron(H1,kron(H1, H2));
elseif s==32
    H = kron(H2,kron(H1, H2));
elseif s==64
    H = kron(H1,kron(H2,kron(H1, H2)));
elseif s==128
    H = kron(H1,kron(H2,kron(H1,kron(H1, H2))));
elseif s==256
    H = kron(H1,kron(H1,kron(H2,kron(H1,kron(H1, H2)))));
elseif s==512
    H = kron(H2,kron(H1,kron(H2,kron(H1,kron(H1, H2)))));
elseif s==1024
    H = kron(H1,kron(H2,kron(H1,kron(H2,kron(H1,kron(H1, H2))))));
elseif s==2048
    H = kron(H2,kron(H2,kron(H1,kron(H2,kron(H1,kron(H1, H2))))));
elseif s==4096
    H = kron(H1,kron(H2,kron(H2,kron(H1,kron(H2,kron(H1,kron(H1, H2)))))));
elseif s==8192
    H = kron(H2,kron(H2,kron(H2,kron(H1,kron(H2,kron(H1,kron(H1, H2)))))));
end

D = zeros(2*nL,1);
t = zeros(2*nL,1);
q = zeros(2*nL,1);
bar=waitbar(0,'Looping over Rademacher vectors...');
for k=1:s
    % generate Hadamard vector of required length on-the-fly
    % (as the nL-by-s matrix H can become too large to store)
    v_k_rep = repmat(H(:,k), ceil(2*nL/s), 1);
    v_k = v_k_rep(1:2*nL);
    
    % evaluate sumandsu using funm_kryl
    t = t + funm_kryl(A_bipartite,v_k,param) .* v_k;
    
    % vector for elementwise division in diagonal estimator
    % (always the one vector for Rademacher vectors)
    q = q + v_k.^2;
    
    % update diagonal estimator
    D = t ./ q;
    waitbar(k/s,bar);
end
close(bar)
toc

% display largest k diagonal entries (joint centralities)
k=50;
[D_sorted, D_sorted_ind] = sort(D, 'descend');
% ['centrality value', 'broadcaster(logical)', 'layer-ID', 'node-ID']
[D_sorted(1:k), D_sorted_ind(1:k)<nL, ceil(mod(D_sorted_ind(1:k),nL)/n), mod(mod(D_sorted_ind(1:k),nL),n)]

%% split diagonal into broadcaster and receiver centralties and save
D_broadcaster = D(1:nL);
D_receiver = D(nL+1:2*nL);
save('results/diagonal_estimation_SC.mat','D_broadcaster','D_receiver');
