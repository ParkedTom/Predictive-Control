function [ Phi, Gamma, Lambda ] = prediction_matrices( A, B, C, N, yk )
%Takes discrete time model matrices A,B,C and calulates the prediction
%matrices Phi, Gamma and Lambda for N samples. takes the index of y, either
%1 or 0.


[n,m] = size(B);

if yk > 0
    Atilda = sparse([A;zeros(n*(N-1),n)]);% calculate Atilda matrix
    Bbar = sparse(kron(eye(N),B));
else
    Atilda = sparse([eye(n);zeros((N-1)*n,n)]);
    Bbar = sparse([zeros(n,(N)*m);kron(eye(N-1),B),zeros((N-1)*n,m)]);
end

Lambda = sparse(kron(eye(N),C));
M = sparse(eye(N*n)) - sparse([zeros(n,n*(N-1)),zeros(n);kron(eye(N-1),A),zeros(n*(N-1),n)]); 


Phi = M\Atilda;
Gamma = M\Bbar;




end

