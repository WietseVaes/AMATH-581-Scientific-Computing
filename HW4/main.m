%% n = 128
clear all
%{
for lambda = 0:0.01:1
    g = @(z) abs(1+lambda/6*(-cos(2*z)+16*cos(z)-15));
    %fplot(g,[pi/2,pi])
    if  g(fminbnd(g,pi/2,pi))<1
        lambdamin = lambda;
    end
end
%}
L = 10;
T = 2;
alpha = 2;
nx = 128;
nt = 501;

x = linspace(-L,L,nx+1); x(end)=[];
Deltax = mean(diff(x));
t = linspace(0,T,nt); Deltat = mean(diff(t));
lambdastar = Deltat*alpha/(Deltax^2);

g1 = @(z)-abs( 1+lambdastar/6*(-cos(2*z)+16*cos(z)-15));
A1 = abs(g1(1));
A2 = abs(g1(fminbnd(g1,-pi,pi)));

e1 = ones(nx,1);
D4 = spdiags([-e1,16*e1,-30*e1,16*e1,-e1],-2:2, nx, nx);
D4(1,nx-1) = -1; D4(1,nx) = 16; D4(2,nx) = -1;
D4(nx-1,1) = -1; D4(nx,1) = 16; D4(nx,2) = -1; 
D4 = D4/12;
A3 = full(D4);

f = @(x) 10*cos(2*pi*x/L)+30*cos(8*pi*x/L);

u1 = zeros(nx,nt);u2 = zeros(nx,nt);u3 = zeros(nx,nt);
u1(:,1) = f(x)';
for index1 = 2:nt
    u1(:,index1) = u1(:,index1-1)+lambdastar*D4*u1(:,index1-1);
    %plot(x,u1(:,index1))
end
A5 = u1(:,end);

g2 = @(z) -abs((1+lambdastar*(cos(z)-1))/(1-lambdastar*(cos(z)-1)));
A6 = abs(g2(fminbnd(g2,-pi,pi)));

e1 = ones(nx,1); 
B = spdiags([-lambdastar*e1/2,e1,-lambdastar*e1/2],-1:1, nx, nx); C = spdiags([lambdastar*e1/2, e1, lambdastar*e1/2],-1:1, nx, nx);
B = B + lambdastar*speye(nx,nx); C = C - lambdastar*speye(nx,nx);
B(1,end) = -lambdastar/2;B(end,1) = -lambdastar/2;
C(1,end) = lambdastar/2;C(end,1) = lambdastar/2;

A7 = full(B); A8 = full(C);
[L,U,P] = lu(B);

u2(:,1) = f(x)';
tic
for index1 = 2:nt
    u2(:,index1) = U\(L\(P*(C*u2(:,index1-1))));
end
toc
A9 = u2(:,end);

u3(:,1) = f(x)';
tic
for index1 = 2:nt
    u3(:,index1) = bicgstab(B,C*u3(:,index1-1));
end
toc
A10 = u3(:,end);

u128 = readmatrix('exact_128.csv');
A11 = norm(u128-A5); A12 = norm(u128-A9);

%% n = 256
clearvars -except A1 A2 A3 A4 A5 A6 A7 A8 A9 A10 A11 A12 lambdastar
L = 10;
T = 2;
alpha = 2;

nx = 256;
x = linspace(-L,L,nx+1); x(end)=[];
Deltax = mean(diff(x));

Deltat = lambdastar/alpha*Deltax^2;

t = 0:Deltat:T; nt = length(t);

e1 = ones(nx,1);
D4 = spdiags([-e1,16*e1,-30*e1,16*e1,-e1],-2:2, nx, nx);
D4(1,nx-1) = -1; D4(1,nx) = 16; D4(2,nx) = -1;
D4(nx-1,1) = -1; D4(nx,1) = 16; D4(nx,2) = -1; 
D4 = D4/12;

f = @(x) 10*cos(2*pi*x/L)+30*cos(8*pi*x/L);

u1 = zeros(nx,nt);u2 = zeros(nx,nt);u3 = zeros(nx,nt);
u1(:,1) = f(x)';
for index1 = 2:nt
    u1(:,index1) = u1(:,index1-1)+lambdastar*D4*u1(:,index1-1);
end

e1 = ones(nx,1); 
B = spdiags([-lambdastar*e1/2,e1,-lambdastar*e1/2],-1:1, nx, nx); C = spdiags([lambdastar*e1/2, e1, lambdastar*e1/2],-1:1, nx, nx);
B = B + lambdastar*speye(nx,nx); C = C - lambdastar*speye(nx,nx);
B(1,end) = -lambdastar/2;B(end,1) = -lambdastar/2;
C(1,end) = lambdastar/2;C(end,1) = lambdastar/2;

[L,U,P] = lu(B);

u2(:,1) = f(x)';
for index1 = 2:nt
    u2(:,index1) = U\(L\(P*(C*u2(:,index1-1))));
end

u256 = readmatrix('exact_256.csv');
A13 = norm(u256-u1(:,end)); A14 = norm(u256-u2(:,end));