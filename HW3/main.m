clear all
%% Problem 1
f = @(x) exp(-(x-5).^2);
L = 10;
T = 10;
Deltax = 0.1;

x = linspace(-L,L-Deltax,200);
N = length(x);

A = sparse(N,N);
for index1 = 1:N
    if index1 ~= 1 && index1 ~= N
        A(index1,index1-1) = -1;A(index1,index1+1) = 1;
    elseif index1 == N
        A(index1,index1-1) = -1;A(index1,1) = 1;
    elseif index1 == 1
        A(index1,index1+1) = 1; A(index1,N) = -1;
    end
end

A = A/(Deltax*2);
A1=full(A);

y0 = f(x);

[t, y1] = ode45(@(t, x) advec1(t, x, A), linspace(0,T,21), y0);
A2 = y1';

xx=x;
[t, y2] = ode45(@(t, x) advec2(t, x, xx), linspace(0,T,21), y0);

A3 = y2';
figure(1);
[XX,TT] = meshgrid(t,x);
contourf(XX,TT,y1','edgecolor','none'); xlabel('t','Interpreter','Latex','FontSize',14);ylabel('x','Interpreter','Latex','FontSize',14); title(sprintf('Solution to the advection equation using \n c=-.5'),...
        'Interpreter','Latex')
figure(2);
contourf(XX,TT,y2','edgecolor','none'); xlabel('t','Interpreter','Latex','FontSize',14);ylabel('x','Interpreter','Latex','FontSize',14); title(sprintf('Solution to the advection equation using \n c=-(1+2sin(5t)-H(x-4))'),...
        'Interpreter','Latex')
hold on; yline(4,'r-','Linewidth',1.5)
%% Problem 2
clearvars -except A1 A2 A3 
m = 64;
n=m*m;
e1 = ones(n,1);
Low1 = repmat([ones(m-1, 1); 0], m, 1); 
Low2 = repmat([1; zeros(m-1, 1)], m, 1); 
Up1 = circshift(Low1, 1);
Up2 = circshift(Low2, m-1);
A = spdiags([e1, e1, Low2, Low1, -4*e1, Up1, Up2, e1, e1], ...
       [-(n-m), -m, -m+1, -1, 0, 1, m-1, m, (n-m)], n, n);

A(1,1)=2;

B = spdiags([e1, -e1, e1, -e1], [-(n-m), -m, m, (n-m)], n, n);

C = spdiags([Low2, -Low1, Up1, -Up2], [-m+1, -1,  1, m-1], n, n);

nu = 0.001;
f = @(x,y)  exp(-2*x.^2-(y.^2/20));
L = 10;
T = 4;
Deltax = 20/64;
Deltat = 0.5;

xx = linspace(-L,L-Deltax,64)'; yy = xx;

A = A/(Deltax^2);
B = B/(2*Deltax);
C = C/(2*Deltax);

A4 = full(A); A5 = full(B); A6 = full(C);

y = repmat(yy,64,1);
x = repmat(xx,1,64)';
x=reshape(x,[64^2,1]);

init = A\f(x,y);
[L,U,P] = lu(A);

tic
[t, ppsi1] = ode45(@(t, x) vorticity1(t, x, A,B,C,nu), 0:Deltat:T, init);
toc

tic
[t, ppsi2] = ode45(@(t, x) vorticity2(t, x, A,B,C,nu,L,U,P), 0:Deltat:T, init);
toc

for index1 = 1:length(t)
    Omega1(index1,:) = A*ppsi1(index1,:)';
    Omega2(index1,:) = A*ppsi2(index1,:)';

end

for index1 = 1:length(t)
     omega1(index1,:,:) = reshape(Omega1(index1,:),[64,64]);
     psi1(index1,:,:) = reshape(ppsi1(index1,:)',[64,64]);
    omega2(index1,:,:) = reshape(Omega2(index1,:),[64,64]);
     psi2(index1,:,:) = reshape(ppsi2(index1,:)',[64,64]);
end

A7 = Omega1; A8 = Omega2;
A9 = omega2;


%{
[XX,YY] = meshgrid(xx,yy);
figure(3);
subplot(1,2,1)
contourf(XX,YY,reshape(psi1(end,:),[64,64]),100,'edgecolor','none');colormap jet
subplot(1,2,2)
contourf(XX,YY,reshape(omega1(end,:),[64,64]),100,'edgecolor','none');colormap jet
figure(4);
subplot(1,2,1)
contourf(XX,YY,reshape(psi2(end,:),[64,64]),100,'edgecolor','none');colormap jet
subplot(1,2,2)
contourf(XX,YY,reshape(omega2(end,:),[64,64]),100,'edgecolor','none');colormap jet
%}

clearvars -except A1 A2 A3 A4 A5 A6 A7 A8 A9

%% Functions
function u_t = advec1(t, x, A)
c = -0.5;
u_t = -c*A*x;
end

function u_t = advec2(t, x, xx)
    N = length(xx);
    A = sparse(N,N);
    Deltax = 0.1;
    for index1 = 1:N
        c = (1-heaviside(xx(index1)-4)+2*sin(5*t));
        if index1 ~= 1 && index1 ~= N
            A(index1,index1-1) = -c;A(index1,index1+1) = c;
        elseif index1 == N
            A(index1,index1-1) = -c;A(index1,1) = c;
        elseif index1 == 1
            A(index1,index1+1) = c; A(index1,N) = -c;
        end
    end
    A = A/(2*Deltax);
    u_t = A*x;
end

function H = heaviside(x)
    H = (x > 0) + 0.5*(x == 0);
end

function ut = vorticity1(t, x, A,B,C,nu)
    b =  nu*A^2*x+(C*x).*(B*A*x)-(B*x).*(C*A*x);
    ut = A\b;
end
function ut = vorticity2(t, x, A,B,C,nu,L,U,P)
    b =  nu*A^2*x+(C*x).*(B*A*x)-(B*x).*(C*A*x);
    ut = U\(L\(P*b));
end