Ac = [-1.2822,0,0.98,0;0,0,1,0;-5.4293,0,-1.8366,0;-128.2,128.2,0,0]; %continuous  time state free response matrix
Bc = [-0.3;0;-17;0]; %continuous time forced response matrix
Cc = [0,1,0,0;0,0,0,1;-128.2,128.2,0,0]; % state-output matrix
Ts = 0.5; % sampling time.


%[A,B] = cont2discrete(Ac,Bc,Ts);
[A,B,C] = cont2discrete(Ac,Bc,Cc,0,Ts);% calculate discrete state space models

n = length(A);

%question 1: return A(4,2)
fprintf('Question 1: %2.1f\n',A(4,2));

%question 2: calculate A(4,2) for Ts = 0.1s
Ts = 0.1;
%[A,B] = cont2discrete(Ac,Bc,Ts);
[A,B,C] = cont2discrete(Ac,Bc,Cc,0,Ts);% calculate discrete state space models
fprintf('Question 2: %2.1f\n',A(4,2));

%question 3: Number of rows of gamma Ts = 0.5s predict 10s into future
Ts = 0.5;
N = 10/0.5;
GammaRowLength = n*N;

fprintf('Question 3: %d\n', GammaRowLength);

%question 4: number of rows of kron(In, C)
N = 16;
[m,n] = size(kron(eye(N),Cc));
fprintf('Question 4: %d\n',m);

%question 5: find gamma(60,2) for Ts = 0.2 and N = 32
Ts = 0.2;
N = 32;
%[A,B] = cont2discrete(Ac,Bc,Ts);
[A,B,C] = cont2discrete(Ac,Bc,Cc,0,Ts);% calculate discrete state space models
[phi,gamma,lambda] = prediction_matrices(A,B,C,N);
fprintf('Question 5: %2.1f\n', full(phi(60,2)));

%question 6: find gamma(192,36) for Ts = 0.1s, N = 64
Ts = 0.1;
N = 64;
%[A,B] = cont2discrete(Ac,Bc,Ts);
[A,B,C] = cont2discrete(Ac,Bc,Cc,0,Ts);% calculate discrete state space models
[phi,gamma,lambda] = prediction_matrices(A,B,C,N);
fprintf('Question 6: %2.1f\n', full(gamma(192,36)));

%question 7: Find 2-norm of Xn after 8s
Ts = 0.5;
N = 16;
X0 = [0;0.1;1;1000];
%Ubar = [0;0;-0.1;-0.1;-0.1;-0.1;0;0;];
Ubar = [0;0;0;0;-0.1;-0.1;-0.1;-0.1;-0.1;-0.1;-0.1;-0.1;0;0;0;0];
%[A,B] = cont2discrete(Ac,Bc,Ts);
[A,B,C] = cont2discrete(Ac,Bc,Cc,0,Ts);% calculate discrete state space models
[phi,gamma,lambda] = prediction_matrices(A,B,C,N);
Xbar = phi*X0 + gamma*Ubar;
n = size(A,2);
norm(Xbar((n*N)-3:(n*N)))




