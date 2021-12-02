% Description: compares the upper bound on resolvent-based subgraph
% centrality computed by Gauss--Radau quadrature applied to the symmetric
% bipartite network representation in dependence of the number of Lanczos 
% iterations to the 'exact' quantity computed with the backslash for the 
% temporal E-Mail EU network and plots the approximation error used in
% [1, Fig. 7].
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

%% assemble symmetric bipartite network representation
A_bipartite=[zeros(nL,nL), A; A', zeros(nL,nL)];

% spy plot of bipartite supra-adjacency matrix
figure(1)
spy(A_bipartite)

% compute smallest and largest eigenvalue of A_bipartite
lambda_min=eigs(A_bipartite,1,'smallestreal');
lambda_max=eigs(A_bipartite,1,'largestreal');

%% compute resovent-based subgraph centrality
% set parameter alpha
alpha=0.5/lambda_max;

% 'exact' result computed with the backslash
SCres=diag((eye(size(A_bipartite))-alpha*A_bipartite)\eye(size(A_bipartite)));

% specify range of iterations
maxit=2:10;
upper_gauss_radau_resolvent=zeros(2*nL,length(maxit));

% create function handle, which realizes the matrix-vector product with
% A_bipartite
MV_prod_bipartite_fh=@(v) MV_prod_bipartite_adjacency(A_layer_adjacencies, interlayer_adjacency_matrix, omegas(1), v);

% loop over range of iterations to approximate SCres with Gauss
% quadrature rules
for j=1:length(maxit)
    bar=waitbar(0,'Looping over nodes');
    tic
    % loop over unit vectors
    for i=1:2*nL
        % set unit vectors as starting vectors
        ei=zeros(2*nL,1); ei(i)=1;
        % Lanczos tridiagonalization using the bipartite function handle
        T_jp1=lanczos_tridiag_fh(MV_prod_bipartite_fh,ei,maxit(j));
        
        % Compute upper Gauss--Radau bound on e_i^T*f(A)*e_i
        upper_gauss_radau_resolvent(i,j)=gauss_radau_resolvent(T_jp1,alpha,lambda_max);
        waitbar(i/(2*nL),bar);
    end
    toc
    close(bar)
end

%% error plot
figure(2)
semilogy(maxit,max(abs(SCres-upper_gauss_radau_resolvent)))

