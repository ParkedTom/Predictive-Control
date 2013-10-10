function [ Phi, Gamma, Lambda ] = prediction_matrices( A, B, C, N )
%Takes discrete time model matrices A,B,C and calulates the prediction
%matrices Phi, Gamma and Lambda.

n = length(A);
Atilda = [A;zeros(n*(N-1),n)];
Bbar = sparse(kron(eye(N),B));
Lambda = sparse(kron(eye(N),C));
M = sparse(eye(N*n)) - sparse([zeros(n,n*(N-1)),zeros(n);kron(eye(N-1),A),zeros(n*(N-1),n)]);
Phi = M\Atilda;
Gamma = M\Bbar;




end

