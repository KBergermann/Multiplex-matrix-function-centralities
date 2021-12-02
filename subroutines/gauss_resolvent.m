function int = gauss_resolvent(T, alpha_resolvent)
% Description: computes lower Gauss quadrature bound on resolvent-based 
% subgraph centrality, i.e., u^T f(A) u for u \in R^n, A \in R^{n \times n}
% symmetric, and f(A) = (I - alpha_resolvent*A)^(-1) the matrix resolvent [1].
% 
% [1] Golub, G. H. & Meurant, G. (2009) Matrices, moments and quadrature 
% with applications. Princeton University Press.
% 
% Input:    T: tridiagonal matrix obtained from the symmetric Lanczos
%               method with the matrix A and starting vector u
%           alpha_resolvent: scalar parameter (between 0 and 1/lambda_max 
%               with lambda_max the largest eigenvalue of A) in the matrix 
%               resolvent function
% Output:   lower bound on u^T f(A) u.
% 
% Kai Bergermann, 2021

j=size(T,1)-1;
[EV,EW] = eig(T(1:j,1:j));
e = zeros(j,1);
e(1) = 1;
int=e'*EV*((eye(j,j)-alpha_resolvent*EW)\eye(j,j))*EV'*e;
end