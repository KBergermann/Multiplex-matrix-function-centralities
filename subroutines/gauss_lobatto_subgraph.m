function int = gauss_lobatto_subgraph(T, beta_subgraph, lambda_min, lambda_max)
% Description: computes upper Gauss--Lobatto quadrature bound on 
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
%           lambda_min: smallest eigenvalue of A
%           lambda_max: largest eigenvalue of A
% Output:   upper bound on u^T f(A) u.
% 
% Kai Bergermann, 2021

j=size(T,1)-1;
e_j=zeros(j,1);
e_j(j)=1;
% prescribe both eigenvalues lambda_min and lambda_max to T
delta=(T(1:j,1:j)-lambda_min*eye(j,j))\e_j;
mu=(T(1:j,1:j)-lambda_max*eye(j,j))\e_j;
T_entries=[1,-delta(j);1,-mu(j)]\[lambda_min;lambda_max];
T(j+1,j+1)=T_entries(1); 
T(j,j+1)=sqrt(T_entries(2)); T(j+1,j)=sqrt(T_entries(2));
e = zeros(j+1,1);
e(1) = 1;
[EV,EW] = eig(T);
int=e'*EV*diag(exp(beta_subgraph*diag(EW)))*EV'*e;
end