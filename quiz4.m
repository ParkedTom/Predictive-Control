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
bla = [Qbar*Sybar - Qbar*lambda* phi*x0;Rbar*Subar];


%contraints
ul = (-15*pi)/180;
uh = (15*pi)/180;
yh = [0.35;410;30];
yl = -yh;

%constraint matrices
Du = [eye(N);-eye(N)];
fu = [kron(ones(N,1),uh);-kron(ones(N,1),ul)];
Dy = [lambda*gamma;-lambda*gamma];
fy = [kron(ones(N,1),yh);kron(-ones(N,1),yl)] - [lambda*phi;-lambda*phi]*x0;
D = [Du;Dy];
fc = [fu;fy];

H = Ala.'*Ala;
f = -Ala.'*bla;

Ubar = quadprog(H,f,D,fc);

x0 = [0;0;0;400];
K = 20/Ts;

Y = C*x0;

Ybar = Y;
Xbar = x0;
Ubar0 = 0;
for k = 1:K
    bla = [Qbar*Sybar - Qbar*lambda*phi*x0;Rbar*Subar];
    f = -Ala.'*bla;
    fy = [kron(ones(N,1),yh);kron(-ones(N,1),yl)] - [lambda*phi;-lambda*phi]*x0;
    fc = [fu;fy];
    Ubar = quadprog(H,f,D,fc);
    x0 = A*x0 + B*Ubar(1);
    %uh - Ubar(1) 
    Ubar(1)
    Y = C*x0;
    Ybar = [Ybar,Y];
    %yh - Y
   
end
Yt = Ybar.'
[a,b] = size(Yt);
n = 0;
for i = 1:a
    if abs(Yt(i,3)- yl(3)) < 0.01
        n = n+1;
    end
end
n

x0 = [0;0.2;0;400];
fy = [kron(ones(N,1),yh);kron(-ones(N,1),yl)] - [lambda*phi;-lambda*phi]*x0;
fc = [fu;fy];
Cinf = [zeros(N,1);1];
Dinf = [D,-ones(80,1);zeros(1,10),-1];
finf = [fc;0];
T = linprog(Cinf,Dinf,finf)

E = 0.0001;
Y = C*x0;

Ybar = Y;
for k = 1:K
    fy = [kron(ones(N,1),yh);kron(-ones(N,1),yl)] - [lambda*phi;-lambda*phi]*x0;
    fc = [fu;fy];
    Cinf = [zeros(N,1);1];
    Dinf = [D,-ones(80,1);zeros(1,10),-1];
    finf = [fc;0];
    T = linprog(Cinf,Dinf,finf);
    bla = [Qbar*Sybar - Qbar*lambda*phi*x0;Rbar*Subar];
    f = -Ala.'*bla;
    fy = [kron(ones(N,1),yh) + ones(3*N,1)*T(11) + E;kron(-ones(N,1),yl) + ones(3*N,1)*T(11) + E] - [lambda*phi;-lambda*phi]*x0;
    fc = [fu;fy];
    Ubar = quadprog(H,f,D,fc);
    x0 = A*x0 + B*Ubar(1);
    Y = C*x0;
    Ybar = [Ybar,Y];
end

