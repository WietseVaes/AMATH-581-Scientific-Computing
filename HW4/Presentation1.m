clear all
%% Parameters + Getting Lambda
L = 10;
T = 2;
alpha = 2;
Allnx = [128,256,512,1024];
nt = 501;

x = linspace(-L,L,Allnx(1)+1); x(end)=[];
Deltax = mean(diff(x));
t = linspace(0,T,nt); Deltat = mean(diff(t));

lambda = Deltat*alpha/(Deltax^2);


%% calculating all approximations
%Calculating aproximations using Method of Lines (1,4) and Crank-Nicolson Method
for index1 = 1:4
    nx = Allnx(index1);
    [u1{index1},u2{index1},DeltaX(index1)] = calcu(nx,lambda);
end

%% calculating erros
%Reading in exact solutions
names = {'exact_128.csv','exact_256.csv','exact_512.csv','exact_1024.csv'};
for index1 = 1:4
    u{index1} = readmatrix(names{index1});
end

%Calculating 2-norm of error
for index1 = 1:4
        err1(index1)= sqrt(trapz(-L:DeltaX(index1):L-DeltaX(index1),(u1{index1}(:,end)-u{index1}).^2));
        err2(index1)= sqrt(trapz(-L:DeltaX(index1):L-DeltaX(index1),(u2{index1}(:,end)-u{index1}).^2));
end
%% plotting
%plot
figure(1); hold on
coeff1 = polyfit(log(DeltaX),log(err1),1);
coeff2 = polyfit(log(DeltaX),log(err2),1);
loglog(DeltaX,exp(coeff1(2))*DeltaX.^coeff1(1),'--','Color',"b")
loglog(DeltaX,exp(coeff2(2))*DeltaX.^coeff2(1),'--','Color',"r")
legend('(1,4)-accurate scheme','Crank-Nicolson Method')
loglog(DeltaX,err1,'*','Color',"b");loglog(DeltaX,err2,'*','Color',"r");
%setting axis + labeling
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(gca, 'YScale', 'log','XScale', 'log');xlim([min(DeltaX)-.001,max(DeltaX)+.01]);xticks(flip(DeltaX));xticklabels(["$$\frac{20}{1024}$$","$$\frac{20}{512}$$","$$\frac{20}{256}$$","$$\frac{20}{128}$$"])
xlabel('$$\Delta x$$','Interpreter','Latex','FontSize',16); ylabel('$$\|u-\tilde{u}\|_2$$','Interpreter','Latex','FontSize',16);
title('Error between heat solution $$u$$ and approximation $$\tilde{u}$$','Interpreter','Latex','FontSize',16)
%legend
lgd = legend('(1,4)-accurate scheme','Crank-Nicolson Method');
lgd.FontSize = 13;

%% Function
function [u1,u2,Deltax] = calcu(nx,lambda)
%Parameters
L = 10;
T = 2;
alpha = 2;
%x-grid
x = linspace(-L,L,nx+1); x(end)=[];
Deltax = mean(diff(x));
%t-grid
Deltat = lambda/alpha*Deltax^2;
t = 0:Deltat:T; nt = length(t);

%making D4
e1 = ones(nx,1);
D4 = spdiags([-e1,16*e1,-30*e1,16*e1,-e1],-2:2, nx, nx);
D4(1,nx-1) = -1; D4(1,nx) = 16; D4(2,nx) = -1;
D4(nx-1,1) = -1; D4(nx,1) = 16; D4(nx,2) = -1; 
D4 = D4/12;

%initial condition
f = @(x) 10*cos(2*pi*x/L)+30*cos(8*pi*x/L);
u1 = zeros(nx,nt);u2 = zeros(nx,nt);
u1(:,1) = f(x)';
for index1 = 2:nt
    u1(:,index1) = u1(:,index1-1)+lambda*D4*u1(:,index1-1); %calculating solution using MOL
end
%Making  B and C for Crank-Nicolson Method
e1 = ones(nx,1); 
B = spdiags([-lambda*e1/2,e1,-lambda*e1/2],-1:1, nx, nx); C = spdiags([lambda*e1/2, e1, lambda*e1/2],-1:1, nx, nx);
B = B + lambda*speye(nx,nx); C = C - lambda*speye(nx,nx);
B(1,end) = -lambda/2;B(end,1) = -lambda/2;
C(1,end) = lambda/2;C(end,1) = lambda/2;
%Calculating the LU-decomposition
[L,U,P] = lu(B);

u2(:,1) = f(x)';
for index1 = 2:nt
    u2(:,index1) = U\(L\(P*(C*u2(:,index1-1)))); %calculating solution using LU
end

end