% Description: compares Gauss quadrature bounds on resolvent-based subgraph
% centrality to the 'exact' quantity computed with the backslash for the
% weighted Scotland Yard network and plots the approximation errors. The
% upper Gauss--Radau bound  is displayed in [1, Fig. 7].
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
A_intra = spalloc(nL,nL,nL);
for t=1:L
  ids = (t-1)*n+(1:n);
  A_intra(ids,ids) = A_layer_adjacencies{t}; %intra-layer block matrix
end

omegas = 10^0;

A = A_intra + omegas(1)*kron(interlayer_adjacency_matrix,speye(n));

% spy plot of supra-adjacency matrix
figure(1)
spy(A)

% compute smallest and largest eigenvalue
lambda_min=eigs(A,1,'smallestreal');
lambda_max=eigs(A,1,'largestreal');

%% compute resolvent-based subgraph centrality
alpha=0.5/lambda_max;

% 'exact' result computed with the backslash
SC_res=diag((eye(nL,nL)-alpha*A)\eye(nL,nL));

% specify range of iterations
maxit=2:10;
lower_gauss_resolvent=zeros(nL,length(maxit));
lower_gauss_radau_resolvent=zeros(nL,length(maxit));
upper_gauss_radau_resolvent=zeros(nL,length(maxit));
upper_gauss_lobatto_resolvent=zeros(nL,length(maxit));

% loop over range of iterations to approximate SCres with Gauss
% quadrature rules
for j=1:length(maxit)
    bar=waitbar(0,'Looping over nodes');

    tic
    for i=1:nL
        % set unit vectors as starting vectors
        ei=zeros(nL,1); ei(i)=1;
        % Lanczos tridiagonalization
        T_jp1=lanczos_tridiag(A,ei,maxit(j));
        
        % Compute lower and upper bounds on e_i^T*f(A)*e_i
        lower_gauss_resolvent(i,j)=gauss_resolvent(T_jp1,alpha);
        lower_gauss_radau_resolvent(i,j)=gauss_radau_resolvent(T_jp1,alpha,lambda_min);
        upper_gauss_radau_resolvent(i,j)=gauss_radau_resolvent(T_jp1,alpha,lambda_max);
        upper_gauss_lobatto_resolvent(i,j)=gauss_lobatto_resolvent(T_jp1,alpha,lambda_min,lambda_max);
        waitbar(i/nL,bar);
    end
    close(bar)
    toc
end

%% error plot
figure(2)
semilogy(maxit,max(abs(SC_res-lower_gauss_resolvent)))
hold on
semilogy(maxit,max(abs(SC_res-lower_gauss_radau_resolvent)))
semilogy(maxit,max(abs(SC_res-upper_gauss_radau_resolvent)))
semilogy(maxit,max(abs(SC_res-upper_gauss_lobatto_resolvent)))
hold off

