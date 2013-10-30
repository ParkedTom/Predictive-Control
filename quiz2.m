A = [0.1,0.2;-1,0.6];
B = [0;1];
C = [0,1;1,0];
x0 = [0;0];
N = 10;
s = [1;0];
Sbar = [s;s;s;s;s;s;s;s;s;s];

[phi,gamma,lambda] = prediction_matrices(A,B,C,N);

Ala = lambda*gamma;
Bla = (Sbar - lambda*phi*x0);

Ubar = Ala.'*Ala\Ala.'*Bla;
Ubar(1)

Ac = [-1.2822,0,0.98,0;0,0,1,0;-5.4293,0,-1.8366,0;-128.2,128.2,0,0]; %continuous  time state free response matrix
Bc = [-0.3;0;-17;0]; %continuous time forced response matrix
Cc = [0,1,0,0;0,0,0,1;-128.2,128.2,0,0]; % state-output matrix
Ts = 0.5; % sampling time.

[A,B,C] = cont2discrete(Ac,Bc,Cc,0,Ts);
Q = eye(3);
R = 1;
x0 = [0;0;0;40];
Sy = [0;0;0];
Su = 0;
N = 10;
n = 4;
Sybar = Sy;
Subar = Su;
for i = 2:N
    Sybar = [Sybar;Sy];
    Subar = [Subar;Su];
end 
   
[phi,gamma,lambda] = prediction_matrices(A,B,C,N-1);
Qbar = kron(eye(N),Q);
Rbar = kron(eye(N),R);

Ala = [Qbar*[zeros(3,10);lambda*gamma,zeros(27,1)];Rbar];
bla = [Qbar*Sybar - Qbar*[zeros(3,4);lambda*phi]*x0;Rbar*Subar];

Ala(23,4)

Ustar = Ala\bla;

(Ustar*180)/pi





