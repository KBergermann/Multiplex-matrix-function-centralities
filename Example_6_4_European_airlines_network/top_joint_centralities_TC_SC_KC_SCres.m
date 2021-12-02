% Description: computes approximations to the top 10 joint centralities of
% total communicability (TC), subgraph centrality (SC), Katz centrality
% (KC), and resolvend-based subgraph centrality (SCres) of the European 
% airlines multiplex network [1,2]. TC and KC are approximated with Krylov 
% subspace methods and for SC and SCres we choose the upper Gauss--Radau 
% bound as approximation. The top 10 node-layer pairs of all measures are 
% returned alongside the airport and airline names corresponding to the 
% involved node-layer pair IDs. The results of KC are displayed and 
% compared with degree and eigenvector centrality [2,3] in [4, Fig. 3].
%  
% [1] Cardillo, A., Gómez-Gardenes, J., Zanin, M., Romance, M., Papo, D.,
% Del Pozo, F. & Boccaletti, S. (2013) Emergence of network features from 
% multiplexity. Sci. Rep., 3(1), 1-6.
% [2] Taylor, D. (2021) Code Release: Supracentrality. Available at 
% https://github.com/taylordr/Supracentrality.
% [3] Taylor, D., Porter, M. A. & Mucha, P. J. (2021) Tunable eigenvector-
% based centralities for multiplex and temporal networks. Multiscale Model.
% Simul., 19(1), 113-147.
% [4] Bergermann, K. & Stoll, M. (2021) Fast computation of matrix 
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

%% compute joint centralities of TC, SC, KC, nad SCres
% set parameters alpha and beta
alpha=0.5/lambda_max;
beta=0.5/lambda_max;

% set the number of Lanczos iterations for Gauss quadrature
maxit=5;

%% compute TC using funm_kryl
tic
param.function = @(A) expm(beta*A);
param.restart_length = 5;
[TC,out_TC] = funm_kryl(A,ones(nL,1),param);
toc

%% compute KC using funm_kryl
param.function = @(A) (eye(size(A))-alpha*A)\eye(size(A));
tic
[KC,out_KC] = funm_kryl(A,ones(nL,1),param);
toc

%% compute SC and SCres
lower_gauss_subgraph=zeros(nL,1);
lower_gauss_radau_subgraph=zeros(nL,1);
upper_gauss_radau_subgraph=zeros(nL,1);
upper_gauss_lobatto_subgraph=zeros(nL,1);
lower_gauss_resolvent=zeros(nL,1);
lower_gauss_radau_resolvent=zeros(nL,1);
upper_gauss_radau_resolvent=zeros(nL,1);
upper_gauss_lobatto_resolvent=zeros(nL,1);

bar=waitbar(0,'Looping over nodes');
tic
% loop over unit vectors
for i=1:nL
    % set unit vectors as starting vectors
    ei=zeros(nL,1); ei(i)=1;
    % Lanczos tridiagonalization
    T_jp1=lanczos_tridiag(A,ei,maxit);
    
    % Compute lower and upper bounds on e_i^T*exp(beta*A)*e_i
    lower_gauss_subgraph(i)=gauss_subgraph(T_jp1,beta);
    lower_gauss_radau_subgraph(i)=gauss_radau_subgraph(T_jp1,beta,lambda_min);
    upper_gauss_radau_subgraph(i)=gauss_radau_subgraph(T_jp1,beta,lambda_max);
    upper_gauss_lobatto_subgraph(i)=gauss_lobatto_subgraph(T_jp1,beta,lambda_min,lambda_max);
    
    % Compute lower and upper bounds on e_i^T*inv(I-alpha*A)*e_i
    lower_gauss_resolvent(i)=gauss_resolvent(T_jp1,alpha);
    lower_gauss_radau_resolvent(i)=gauss_radau_resolvent(T_jp1,alpha,lambda_min);
    upper_gauss_radau_resolvent(i)=gauss_radau_resolvent(T_jp1,alpha,lambda_max);
    upper_gauss_lobatto_resolvent(i)=gauss_lobatto_resolvent(T_jp1,alpha,lambda_min,lambda_max);

    waitbar(i/nL,bar);
end
close(bar)
toc

% choose upper Gauss--Radau bound as approximation
SC=upper_gauss_radau_subgraph;
SCres=upper_gauss_radau_resolvent;

%% top 10 joint centrality rankings
[TC_sorted,TC_ind]=sort(TC,'descend');
TC_ranking=[mod(TC_ind(1:10),n),ceil(TC_ind(1:10)/n),TC(TC_ind(1:10))];

[SC_sorted,SC_ind]=sort(SC,'descend');
SC_ranking=[mod(SC_ind(1:10),n),ceil(SC_ind(1:10)/n),SC(SC_ind(1:10))];

[KC_sorted,KC_ind]=sort(KC,'descend');
KC_ranking=[mod(KC_ind(1:10),n),ceil(KC_ind(1:10)/n),KC(KC_ind(1:10))];

[SCres_sorted,SCres_ind]=sort(SCres,'descend');
SCres_ranking=[mod(SCres_ind(1:10),n),ceil(SCres_ind(1:10)/n),SCres(SCres_ind(1:10))];

%% print top 10 joint centrality rankings
fprintf('Top 10 joint centralities:\n')
fprintf('|         TC         |        SC          |         KC         |       SCres        |\n')
fprintf('|-----------------------------------------------------------------------------------|\n')
for i=1:10
    fprintf('|  (%3d,%2d)  %.4f  |  ',int16(TC_ranking(i,1)),int8(TC_ranking(i,2)),TC_ranking(i,3))
    fprintf('(%3d,%2d)  %.4f  |  ',int16(SC_ranking(i,1)),int8(SC_ranking(i,2)),SC_ranking(i,3))
    fprintf('(%3d,%2d)  %.4f  |  ',int16(KC_ranking(i,1)),int8(KC_ranking(i,2)),KC_ranking(i,3))
    fprintf('(%3d,%2d)  %.4f  |\n',int16(SCres_ranking(i,1)),int8(SCres_ranking(i,2)),SCres_ranking(i,3))
end

%% get airport and airline name corresponding to IDs in table
layer_id_list=ceil(KC_ind(1:10)/n); layer_id_list(end+1)=10;
node_id_list=mod(KC_ind(1:10),n); node_id_list(end+1)=83;
names{1,:}=airport_names(node_id_list);
names{2,:}=layer_ids(layer_id_list)';
[num2cell(node_id_list),num2cell(layer_id_list),names{1,:},names{2,:}]
