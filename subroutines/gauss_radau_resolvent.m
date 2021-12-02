function int = gauss_radau_resolvent(T, alpha_resolvent, lambda_bound)
% Description: computes lower or upper Gauss--Radau quadrature bound on 
% resolvent-based subgraph centrality, i.e., u^T f(A) u for u \in R^n, 
% A \in R^{n \times n} symmetric, and f(A) = (I - alpha_resolvent*A)^(-1) 
% the matrix resolvent [1].
% 
% [1] Golub, G. H. & Meurant, G. (2009) Matrices, moments and quadrature 
% with applications. Princeton University Press.
% 
% Input:    T: tridiagonal matrix obtained from the symmetric Lanczos
%               method with the matrix A and starting vector u
%           alpha_resolvent: scalar parameter (between 0 and 1/lambda_max 
%               with lambda_max the largest eigenvalue of A) in the matrix 
%               resolvent function
%           lambda_bound: smallest or largest eigenvalue of A 
%               (for f(A) = (I - alpha_resolvent*A)^(-1)  lambda_min yields
%               a lower and lambda_max an upper bound)
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
int=e'*EV*((eye(j+1,j+1)-alpha_resolvent*EW)\eye(j+1,j+1))*EV'*e;
end
