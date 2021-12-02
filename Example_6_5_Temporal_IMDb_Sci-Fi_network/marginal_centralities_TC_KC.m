% Description: computes approximations of broadcaster and receiver total 
% communicability and Katz centrality via Krylov subspace methods for the 
% temporal IMDb Sci-Fi network. The estimates are used to generate the 
% marginal node rankings in [1, Tab. 4].
%
% [1] Bergermann, K. & Stoll, M. (2021) Fast computation of matrix 
% function-based centrality measures for layer-coupled multiplex networks. 
% https://arxiv.org/abs/2104.14368v3
% 
% Kai Bergermann, 2021

%% load subroutines
addpath('../subroutines')
addpath('../subroutines/funm_kryl')

%% load intra-layer adjacency matrices and principal name table
load data/sci_fi_seasons_intra_layer_adjacency.mat
names = readtable('data/names_sci_fi.csv');

% extract network size
L = size(layer_id_table,1);
n = size(node_id_table,1);
nL = n*L;

%% assemble supra-adjacency matrix
% temporal inter-layer coupling with weights exp(-deltaT)
deltaT = layer_id_table(2:end,2)-layer_id_table(1:end-1,2);
inter_layer_adjacency_matrix = spdiags([zeros(L,1),[exp(-deltaT);1]],0:1,L,L);

omega=10^(1);

A = A_intra + omega * kron(inter_layer_adjacency_matrix,speye(n));

%% estimating the extremal eigevalues via an initial lanczos procedure
beta=1e-02;

% set funm_kryl parameters
param.function = @(A) expm(beta*A);        % other choices: 'expBA', 'expCF', ...
param.restart_length = 5;
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
param.waitbar = 1;            % show waitbar 
param.reorth_number = 0;      % #reorthogonalizations
param = param_init(param);    % check and correct param structure

[~,out] = funm_kryl(A,ones(nL,1),param);
lambda=eig(out.H_full);

lambda_min=min(lambda);
lambda_max=max(lambda);

% set parameters alpha and beta
alpha=0.9/lambda_max;
beta=5/lambda_max;

%% compute broadcaster centralities

%% TC
tic
param.function = @(A) expm(beta*A);
param.restart_length = 5;

[TC_broadcast,out_TC_broadcast] = funm_kryl(A,ones(nL,1),param);
toc

%% KC
param.function = @(A) (eye(size(A))-alpha*A)\eye(size(A));
tic
[KC_broadcast,out_KC_broadcast] = funm_kryl(A,ones(nL,1),param);
toc

%% compute top k marginal layer and node centralities
% set number of marginal centralities to display
k=300;

% TC, broadcaster
% compute MLC and MNC
TC_broadcast_MLC = sum(reshape(TC_broadcast,n,L),1);
TC_broadcast_MNC = sum(reshape(TC_broadcast,n,L),2);

% sort TC MNCs and identify corresponding principal names
[~,TC_broadcast_MNC_ind]=sort(TC_broadcast_MNC,'descend');
TC_broadcast_MNC_ranking_names = strings(k,1);
for i=1:k
    TC_broadcast_MNC_ranking_names(i) = string(table2array(names(find(table2array(names(:,3))==node_id_table(mod(TC_broadcast_MNC_ind(i),n),2)),4)));
end
TC_broadcast_MNC_ranking=[TC_broadcast_MNC_ind(1:k),TC_broadcast_MNC_ranking_names,TC_broadcast_MNC(TC_broadcast_MNC_ind(1:k))]

% KC, broadcaster
% compute MLC and MNC
KC_broadcast_MLC = sum(reshape(KC_broadcast,n,L),1);
KC_broadcast_MNC = sum(reshape(KC_broadcast,n,L),2);

% sort KC MNCs and identify corresponding principal names
[~,KC_broadcast_MNC_ind]=sort(KC_broadcast_MNC,'descend');
KC_broadcast_MNC_ranking_names = strings(k,1);
for i=1:k
    KC_broadcast_MNC_ranking_names(i) = string(table2array(names(find(table2array(names(:,3))==node_id_table(mod(KC_broadcast_MNC_ind(i),n),2)),4)));
end
KC_broadcast_MNC_ranking=[KC_broadcast_MNC_ind(1:k),KC_broadcast_MNC_ranking_names,KC_broadcast_MNC(KC_broadcast_MNC_ind(1:k))]

%% compute receiver centralities
%% TC
tic
param.function = @(A) expm(beta*A);
param.restart_length = 5;

[TC_receive,out_TC_receive] = funm_kryl(A',ones(nL,1),param);
toc

%% KC
param.function = @(A) (eye(size(A))-alpha*A)\eye(size(A));
tic
[KC_receive,out_KC_receive] = funm_kryl(A',ones(nL,1),param);
toc

%% compute top k marginal layer and node centralities
% set number of marginal centralities to display
k=300;

% TC, receiver
% compute MLC and MNC
TC_receive_MLC = sum(reshape(TC_receive,n,L),1);
TC_receive_MNC = sum(reshape(TC_receive,n,L),2);

% sort TC MNCs and identify corresponding principal names
[~,TC_receive_MNC_ind]=sort(TC_receive_MNC,'descend');
TC_receive_MNC_ranking_names = strings(k,1);
for i=1:k
    TC_receive_MNC_ranking_names(i) = string(table2array(names(find(table2array(names(:,3))==node_id_table(mod(TC_receive_MNC_ind(i),n),2)),4)));
end
TC_receive_MNC_ranking=[TC_receive_MNC_ind(1:k),TC_receive_MNC_ranking_names,TC_receive_MNC(TC_receive_MNC_ind(1:k))]

% KC, receiver
% compute MLC and MNC
KC_receive_MLC = sum(reshape(KC_receive,n,L),1);
KC_receive_MNC = sum(reshape(KC_receive,n,L),2);

% sort KC MNCs and identify corresponding principal names
[~,KC_receive_MNC_ind]=sort(KC_receive_MNC,'descend');
KC_receive_MNC_ranking_names = strings(k,1);
for i=1:k
    KC_receive_MNC_ranking_names(i) = string(table2array(names(find(table2array(names(:,3))==node_id_table(mod(KC_receive_MNC_ind(i),n),2)),4)));
end
KC_receive_MNC_ranking=[KC_receive_MNC_ind(1:k),KC_receive_MNC_ranking_names,KC_receive_MNC(KC_receive_MNC_ind(1:k))]
