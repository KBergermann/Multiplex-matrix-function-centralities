% Description: compares the approximation of total communicability obtained
% with Krylov subspace methods to the 'exact' quantity computed with expm
% for the temporal E-Mail EU network and plots the approximation
% error for broadcaster and receiver centralities. The maximum error over
% both quantities is plotted in [1, Fig. 7].
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
  A_intra(ids,ids) = A_layer_adjacencies{l};
end

omegas = 10^[0];

A = A_intra + omegas(1) * kron(interlayer_adjacency_matrix,speye(n));

% spy plot of supra-adjacency matrix
figure(1)
spy(A)

% compute largest eigenvalue
lambda_max=eigs(A,1);

%% compute total communicability
% set parameter beta
beta=0.5/lambda_max;

% 'exact' result computed with expm
TC_broadcast=expm(beta*A)*ones(nL,1);
TC_receive=expm(beta*A')*ones(nL,1);

% specify range of iterations
maxit=2:10;
TC_funm_kryl_broadcast=zeros(nL,length(maxit));
TC_funm_kryl_receive=zeros(nL,length(maxit));

% loop over range of iterations to approximate TC with Krylov subspace
% methods
for i=1:length(maxit)
    tic
    % set funm_kryl parameters
    param.function = @(A) expm(beta*A);
    param.restart_length = maxit(i);
    param.max_restarts = 1;
    param.hermitian = 0;          % set 0 if A is not Hermitian
    param.V_full = 0;             % set 1 if you need Krylov basis
    param.H_full = 1;             % if using rational functions you can set this 0
    param.exact = [];          % if not known set to []
    param.bound = 1;              % returns upper and lower bounds (after some cycles)
    param.stopping_accuracy = 1e-12;  % stopping accuracy
    param.inner_product = @inner_product;
    param.thick = [];             % thick-restart function  
    param.min_decay = 0.95;       % we desire linear error reduction of rate < .95 
    param.waitbar = 0;            % show waitbar 
    param.reorth_number = 0;      % #reorthogonalizations
    param = param_init(param);    % check and correct param structure
    
    % approximate TC according to [1, Sec. 5.1]
    [TC_funm_kryl_broadcast(:,i),out1]=funm_kryl(A,ones(nL,1),param);
    [TC_funm_kryl_receive(:,i),out2]=funm_kryl(A',ones(nL,1),param);
    toc
end

%% error plot
figure(2)
semilogy(maxit,max(abs(TC_broadcast-TC_funm_kryl_broadcast)))
hold on
semilogy(maxit,max(abs(TC_receive-TC_funm_kryl_receive)))
hold off

