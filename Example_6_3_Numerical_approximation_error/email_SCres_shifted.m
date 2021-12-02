% Description: compares the approximation of resolvent-based subgraph 
% centrality obtained with the shifted Arnoldi approach [1, Sec. 5.2.2]
% to the 'exact' quantity computed with the backslash for the temporal 
% E-Mail EU network and plots the approximation error used in [1, Fig. 7].
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
load data/email_adjacencies_15d.mat
A_layer_adjacencies = A_intra;

% extract network size
n=size(A_layer_adjacencies{1},1);
L=size(A_layer_adjacencies,1);
nL=n*L;

%% choose inter-layer coupling
interlayer_type=3;
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
for l=1:L
  ids = (l-1)*n + (1:n);
  A_intra(ids,ids) = A_layer_adjacencies{l}; %intra-layer block matrix
end

omegas = 10^[0];

A = A_intra + omegas(1) * kron(interlayer_adjacency_matrix,speye(n));

% spy plot of supra-adjacency matrix
figure(1)
spy(A)

% compute largest eigenvalue
lambda_max=eigs(A,1);

%% compute resolvent-based subgraph centrality
% set parameter alpha
alpha=0.5/lambda_max;

% 'exact' result computed with the backslash
SCres=diag((eye(size(A))-alpha*A)\eye(size(A)));

% specify range of iterations
maxit=2:10;
SCres_shifted=zeros(nL,length(maxit));

% loop over range of iterations to approximate SCres with Krylov subspace
% methods
for j=1:length(maxit)
    %% SC part two: shifted term
    tic
    % set funm_kryl parameters    
    param.function = @(A) (eye(size(A))-alpha*A)\eye(size(A));
    param.restart_length = maxit(j);
    param.max_restarts = 1;
    param.hermitian = 0;          % set 0 if A is not Hermitian
    param.V_full = 0;             % set 1 if you need Krylov basis
    param.H_full = 1;             % if using rational functions you can set this 0
    param.exact = [];             % if not known set to []
    param.bound = 0;              % returns upper and lower bounds (after some cycles)
    param.stopping_accuracy = 1e-12;  % stopping accuracy
    param.inner_product = @inner_product;
    param.thick = [];             % thick-restart function  
    param.min_decay = 0.95;       % we desire linear error reduction of rate < .95 
    param.waitbar = 0;            % show waitbar 
    param.reorth_number = 0;      % #reorthogonalizations
    param = param_init(param);    % check and correct param structure

    % approximate KC according to [1, Sec. 5.1]
    [SCres_part_two,~] = funm_kryl(A,ones(nL,1),param);
    toc

    %% SC part one: e_i^T f(A) [e_i+1]
    SCres_part_one=zeros(nL,1);

    tic
    % loop over unit vectors
    for i=1:nL
        ei=zeros(nL,1); ei(i)=1;
        shift_vector=ones(nL,1);

        % approximate shifted parts according to [1, Sec. 5.1]
        [SCres_part_one_temp,out2] = funm_kryl(A,ei+shift_vector,param);
        SCres_part_one(i)=SCres_part_one_temp(i);
    end
    toc

    % total approximation as described in [1, Sec. 5.2.2]
    SCres_shifted(:,j)=SCres_part_one-SCres_part_two;
end

%% error plot
figure(2)
semilogy(maxit,max(abs(SCres-SCres_shifted)))

