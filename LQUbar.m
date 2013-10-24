function [ Ubar ] = LQUbar( Q, R, Sy, Su, phi, gamma, lambda, x0, N )
%This function calculates Qbar, Rbar, Sybar, Subar and calculates the LQ
%input vector

Qbar = sparse(kron(eye(N),Q));
Rbar = sparse(kron(eye(N),R));
Sybar = Sy;
Subar = Su;
for i = 2:N
    Sybar = [Sybar;Sy];
    Subar = [Subar;Su];
end

Ala = [Qbar*lambda*gamma;Rbar];
bla = [Qbar*Sybar - Qbar*lambda*phi*x0;Rbar*Subar];

Ubar = Ala.'*Ala\Ala.'*bla;

end

