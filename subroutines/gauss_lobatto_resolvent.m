function int = gauss_lobatto_resolvent(T, alpha_resolvent, lambda_min, lambda_max)
% Description: computes upper Gauss--Lobatto quadrature bound on
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
int=e'*EV*((eye(j+1,j+1)-alpha_resolvent*EW)\eye(j+1,j+1))*EV'*e;
end
