% Description: computes deterministic estimates of the trace and the
% diagonal of f(A)=exp(beta*A) using Hadamard vetors in order to obtain
% approximations to the Estrada index and joint subgraph centralities of 
% the European airlines multiplex network [1,2]. The results for different 
% values of s (number of Hadamard vectors) are compared with the explicitly
% computed quantities using expm and are displayed in [3, Fig. 8].
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

%% load intra-layer adjacency matrices
load data/multiplex_airlines_GC.mat

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

%% subgraph centrality explicitly
% set parameter beta
beta = 5/max(lambda);

% 'exact' result computed with expm
expm_betaA = expm(beta*A);
diag_expA = diag(expm_betaA);

% update parameter beta in funm_kryl
param.function = @(A) expm(beta*A);        % other choices: 'expBA', 'expCF', ...

%% subgraph centrality diagonal estimator
% using deterministic Hadamard vectors

tic
% build two small Hadamard matrices from which to compute larger ones via
% Kronecker products
H1 = [1, 1; 1, -1];
H2 = [1, 1, 1, 1; 1, -1, 1, -1; 1, 1, -1, -1; 1, -1, -1, 1];

% set number of Hadamard vectors
s=16;

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
end

% stack H upon itself as often as it requires to have nL rows in H_sy
H_sy = repmat(H, ceil(nL/s), 1);

D = zeros(nL,1);
t = zeros(nL,1);
q = zeros(nL,1);
% loop over Hadamard vectors
for k=1:s
    % generate Hadamard vector of required length
    v_k = H_sy(1:nL,k);
    
    % evaluate sumands using funm_kryl
    t = t + funm_kryl(A,v_k,param) .* v_k;
    
    % vector for elementwise division in diagonal estimator
    % (always the one vector for Rademacher vectors)
    q = q + v_k.^2;
    
    % update diagonal estimator
    D = t ./ q;
end

% compute and output l2-error of diagonal estimation and relative error
% of the trace for 's'-many Hadamard vectors
error_2norm = norm(D-diag_expA,2)/norm(diag_expA,2)
rel_error_EE = abs(sum(D) - sum(diag_expA))/abs(sum(diag_expA))
toc
