% Description: computes approximations to marginal layer and marginal node
% Katz centralities of the European airlines multiplex network [1,2] with
% the coupling parameter \omega varying between 10^(-2) and 10^3. The 
% results are plotted in a loglog-plot and displayed in [3, Fig. 9].
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

%% assemble intra-layer adjacency matrix
A_intra = spalloc(nL,nL,nL);
for t=1:L
  ids = (t-1)*n + (1:n);
  A_intra(ids,ids) = net.A{t}; %intra-layer block matrix
end

%% set range of omegas
omegas = 10.^[-2:0.25:3]; 

%% loop over different values of omega
KC=zeros(nL,length(omegas));
KC_matrix=zeros(n,L,length(omegas));
KC_MLC=zeros(length(omegas),L);
KC_MNC=zeros(n,length(omegas));

for i=1:length(omegas)
    % assemble supra-adjacency matrix
    A = A_intra + omegas(i) * kron(interlayer_adjacency_matrix,speye(n));

    % estimating the extremal eigevalues via an initial lanczos procedure
    alpha=1e-02;

    param.function = @(A) (eye(size(A))-alpha*A)\eye(size(A));
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

    % compute Katz centrality using funm_kryl
    alpha=.999/lambda_max;

    tic
    param.function = @(A) (eye(size(A))-alpha*A)\eye(size(A));
    param.restart_length = 5;
    
    [KC(:,i),out] = funm_kryl(A,ones(nL,1),param);
    toc
    
    % compute marginal layer and node centralities
    KC_matrix(:,:,i)=reshape(KC(:,i),n,L);
    KC_MLC(i,:)=sum(KC_matrix(:,:,i));
    KC_MNC(:,i)=sum(KC_matrix(:,:,i),2);

end

%% plot marginal centrality values
figure(1)
% marginal layer centralities
subplot(121)
loglog(omegas,KC_MLC);
title('Marginal layer centralities')
xlabel('$\omega$','interpreter','latex')
ylabel('MLC(l)')

% marginal node centralities
subplot(122)
loglog(omegas,KC_MNC);
title('Marginal Node Centralities')
xlabel('$\omega$','interpreter','latex')
ylabel('MNC(i)')
