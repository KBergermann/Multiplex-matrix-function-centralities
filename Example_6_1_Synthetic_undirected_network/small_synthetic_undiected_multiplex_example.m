% Description: This script creates a small synthetic multiplex network with
% n=4 nodes and L=3 layers as well as its aggregated single-layer version
% and computes degree as well as matrix function-based centrality measures
% for both networks. The results from [1, Fig. 5] are reproduced with the
% parameters alpha_factor = 0.9 and beta_factor = 3 and the different
% factors \tilde{A}_{12} / \tilde{A}_{13} are obtained by the appropriate
% quantity a12/a13 where the smaller value is set to 1.
% 
% [1] Bergermann, K. & Stoll, M. (2021) Fast computation of matrix 
% function-based centrality measures for layer-coupled multiplex networks. 
% https://arxiv.org/abs/2104.14368v3
% 
% Kai Bergermann, 2021


% specify network size
n=4; 
L=3;

% build intra-layer adjacencies
A1 = zeros(n,n);
A1(1,:) = 1;
A1(:,1) = 1;
A1(1,1) = 0;
A2 = A1;
A2(1,4) = 0; A2(4,1) = 0;
A3 = A1;
A3(1,3) = 0; A3(3,1) = 0;

% set range of parameters alpha and beta via
% alpha = alpha_factor/lambda_max and beta = beta_factor/lambda_max
alpha_factor = 0.9; % [0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99, 0.999];
beta_factor = 3; % [0.01, 0.1, 0.25, 0.5, 1, 3, 5, 20];

% set inter-layer weights \tilde{A}_{12} and \tilde{A}_{13}
a12=10; %sqrt(10);
a13=1;

% loops over specified range of parameters
for i=1:size(alpha_factor,2)
    
    % build supra-adjacency matrix of multiplex network
    interlayer_adjacency_matrix = [0,a12,a13; a12,0,0; a13,0,0];
    omega = 1;
    A = blkdiag(A1,A2,A3) + kron(interlayer_adjacency_matrix, eye(n));

    % compute matrix function-based centralities of multiplex network
    degree_MNC = sum(reshape(sum(A,1),[n,L]),2);
    [V, D] = eig(A);
    lmax=max(max(D));
    alpha=alpha_factor(i)/lmax;
    SCres_MNC = sum(reshape(diag(inv(eye(size(A))-alpha*A)),[n,L]),2);
    KC_MNC = sum(reshape(sum(inv(eye(size(A))-alpha*A),2),[n,L]),2);
    beta= beta_factor(i)/lmax;
    SC = sum(reshape(diag(expm(beta*A)),[n,L]),2);
    TC = sum(reshape(sum(expm(beta*A),2),[n,L]),2);

    % build adjacency matrix of aggregated network
    A_agg = diag(ones(n,1)*2*(a12+a13)) + A1 + A2 + A3;
    
    % compute matrix function-based centralities of aggregated network
    degree_MNC_agg = sum(A_agg,1)';
    [V_agg, D_agg] = eig(A_agg);
    lmax_agg=max(max(D_agg));
    alpha_agg = alpha_factor(i)/lmax_agg;
    SCres_agg = diag(inv(eye(size(A_agg))-alpha_agg*A_agg));
    KC_agg = sum(inv(eye(size(A_agg))-alpha_agg*A_agg),2);
    beta_agg = beta_factor(i)/lmax_agg;
    SC_agg = diag(expm(beta_agg*A_agg));
    TC_agg = sum(expm(beta_agg*A_agg),2);

    % print results
    fprintf('\n Marginal node centrality values for alpha=%.2f/lambda_max, beta=%.2f/lambda_max:', alpha_factor, beta_factor)
    fprintf('\n\n  [degree,    SCres,   KC,       SC,       TC     ]\n')
    multiplex = [degree_MNC, SCres_MNC, KC_MNC, SC, TC]
    fprintf('\n\n  [degree,    SCres,   KC,       SC,       TC     ]\n')
    aggregated = [degree_MNC_agg, SCres_agg, KC_agg, SC_agg, TC_agg]
end
