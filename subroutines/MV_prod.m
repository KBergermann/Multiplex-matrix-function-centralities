function Ab = MV_prod(A_blkdiag, A_tilde, omega, b)
% Description: Function that performs the matrix-vector product A*b
% where A is given by:
% A = blkdiag(A_blkdiag{1}, ... , A_blkdiag{L}) + omega * kron(A_tilde, ones(n,n))
% 
% Input:    A_blkdiag: cell(L,1) with nxn intralayer adjacencies
%           A_tilde: LxL interlayer weight matrix
%           omega: scalar coupling parameter
%           b: vector in 2nL
% Output:   matrix vector product A*b.
% 
% Martin Stoll, Kai Bergermann, 2021

n=size(A_blkdiag{1},1);
L=size(A_tilde,1);

% interlayer multiplication
for i=1:L
    A_blkdiagb((i-1)*n+1:(i)*n,1) = A_blkdiag{i}*b((i-1)*n+1:(i)*n,1);
end
% intralayer multiplication
A_tildeb=reshape(reshape(b,n,L)*(omega*A_tilde)',n*L,1);

Ab=A_blkdiagb+A_tildeb;
end
