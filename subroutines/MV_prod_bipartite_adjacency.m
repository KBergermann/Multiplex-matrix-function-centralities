function Ab = MV_prod_bipartite_adjacency(A_blkdiag, A_tilde, omega, b)
% Description: Function handle that performs the matrix-vector product
% [0, A; A', 0]*b where A is given by:
% A = blkdiag(A_blkdiag{1}, ... , A_blkdiag{L}) + omega * kron(A_tilde, ones(n,n))
% 
% Input:    A_blkdiag: cell(L,1) with nxn intralayer adjacencies
%           A_tilde: LxL interlayer weight matrix
%           b: vector in 2nL
% Output:   matrix vector product [0, A; A', 0]*b.
% 
% Kai Bergermann, 2021

n=size(A_blkdiag{1},1);
L=size(A_tilde,1);

Ab=zeros(2*n*L,1);
Ab(1:n*L)=MV_prod(A_blkdiag, A_tilde, omega, b((n*L+1):end));
Ab(n*L+1:end)=MV_prod_transposed(A_blkdiag, A_tilde, omega, b(1:n*L));
end