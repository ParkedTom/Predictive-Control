Ac = [-1.2822,0,0.98,0;0,0,1,0;-5.4293,0,-1.8366,0;-128.2,128.2,0,0]; %continuous  time state free response matrix
Bc = [-0.3;0;-17;0]; %continuous time forced response matrix
Cc = [0,1,0,0;0,0,0,1;-128.2,128.2,0,0]; % state-output matrix
Ts = 0.5; % sampling time.
N = 10;
[A,B,C] = cont2discrete(Ac,Bc,Cc,0,Ts);

Q = eye(3);
R = 1;
x0 = [0;0;0;400];
Sy = [0;0;0];
Su = 0;

Sybar = kron(ones(N,1),Sy);
Subar = kron(ones(N,1),Su);
[phi,gamma,lambda] = prediction_matrices(A,B,C,N,0);              

Qbar = kron(eye(N),Q);
Rbar = kron(eye(N),R);

Ala = [Qbar*lambda*gamma;Rbar];
bla = [Qbar*Sybar - Qbar*lambda*phi*x0;Rbar*Subar];



ul = (-15*pi)/180;
uh = (15*pi)/180;

%constraints
Du = [eye(N);-eye(N)];
fu = [kron(ones(N,1),uh);-kron(ones(N,1),ul)];

H = Ala.'*Ala;
f = -Ala.'*bla;

Ubar = quadprog(H,f,Du,fu);
(Ubar*180)/pi
x0 = [0;0;0;400];
K = 20/Ts;

Y = C*x0;
if Y(2) < 0.5
    k = 1;
else
    for k = 1:K
        bla = [Qbar*Sybar - Qbar*lambda*phi*x0;Rbar*Subar];
        f = -Ala.'*bla;
        Ubar = quadprog(H,f,Du,fu);
        x0 = A*x0 + B*Ubar(1);
        Y = C*x0;
        if Y(2) < 0.5
            k+1
            break;
        end
    end
end



