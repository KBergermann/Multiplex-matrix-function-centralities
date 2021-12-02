function int = gauss_subgraph(T, beta_subgraph)
% Description: computes lower Gauss quadrature bound on subgraph centrality,
% i.e., u^T f(A) u for u \in R^n, A \in R^{n \times n} symmetric,
% and f(A) = exp(beta_subgraph*A) the matrix exponential [1].
% 
% [1] Golub, G. H. & Meurant, G. (2009) Matrices, moments and quadrature 
% with applications. Princeton University Press.
% 
% Input:    T: tridiagonal matrix obtained from the symmetric Lanczos
%               method with the matrix A and starting vector u
%           beta_subgraph: scalar parameter (inverse temperature) in the
%               matrix exponential
% Output:   lower bound on u^T f(A) u.
% 
% Kai Bergermann, 2021

j=size(T,1)-1;
[EV,EW] = eig(T(1:j,1:j));
e = zeros(j,1);
e(1) = 1;
int=e'*EV*diag(exp(beta_subgraph*diag(EW)))*EV'*e;
end