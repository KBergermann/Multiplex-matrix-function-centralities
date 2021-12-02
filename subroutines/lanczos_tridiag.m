function T = lanczos_tridiag(A,u,maxit)
% Description: Implementation of the symmetric Lanczos method [1]
% 
% [1] Golub, G. H. & Van Loan, C. F. (2013) Matrix computations, volume 3. JHU press.
% 
% Input:    A: symmetric nxn matrix
%           u: starting vector of length n
%           maxit: maximal number of Lanczos iterations (the method 
%               terminates before if stopping criterion is met)
% Output:   tridiagonal matrix T.
% 
% Martin Stoll, Kai Bergermann, 2008-2021

n = length(A(:,1));
u = u/norm(u);
U(:,1) = u;
alpha = [];
beta = [];
T = [];
for j = 1:maxit
    if j == 1
        U(:,j+1) = A*U(:,j);
    else
        U(:,j+1) = A*U(:,j)-beta(j-1)*U(:,j-1);
    end
    alpha(j) = U(:,j+1)'*U(:,j);
    U(:,j+1) = U(:,j+1)-alpha(j)*U(:,j);
    beta(j) = norm(U(:,j+1));
    U(:,j+1) = U(:,j+1)/beta(j);
    T(j,j) = alpha(j);
    T(j,j+1) = beta(j);
    T(j+1,j) = beta(j);
    if abs(beta(j))<1e-12
        break
    end
end
end