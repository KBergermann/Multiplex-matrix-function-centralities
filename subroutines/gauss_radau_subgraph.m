function int = gauss_radau_subgraph(T, beta_subgraph, lambda_bound)
% Description: computes lower or upper Gauss--Radau quadrature bound on 
% subgraph centrality, i.e., u^T f(A) u for u \in R^n, A \in R^{n \times n}
% symmetric, and f(A) = exp(beta_subgraph*A) the matrix resolvent [1].
% 
% [1] Golub, G. H. & Meurant, G. (2009) Matrices, moments and quadrature 
% with applications. Princeton University Press.
% 
% Input:    T: tridiagonal matrix obtained from the symmetric Lanczos
%               method with the matrix A and starting vector u
%           beta_subgraph: scalar parameter (inverse temperature) in the
%               matrix exponential
%           lambda_bound: smallest or largest eigenvalue of A 
%               (for f(A) = exp(beta_subgraph*A) lambda_min yields a lower
%               and lambda_max an upper bound)
% Output:   lower or upper bound on u^T f(A) u.
% 
% Kai Bergermann, 2021

j=size(T,1)-1;
rhs=zeros(j,1);
rhs(j)=T(j,j+1)^2;
e = zeros(j+1,1);
e(1) = 1;
% prescribe the eigenvalue lambda_bound to T
delta_j_lower=(T(1:j,1:j)-lambda_bound*eye(j,j))\rhs;
T(j+1,j+1)=lambda_bound+delta_j_lower(j);
[EV,EW] = eig(T);
int=e'*EV*diag(exp(beta_subgraph*diag(EW)))*EV'*e;
end