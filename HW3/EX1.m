%% Problem 1
clear all
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
        'Interpreter','Latex'); colorbar();
figure(2);
contourf(XX,TT,y2','edgecolor','none'); xlabel('t','Interpreter','Latex','FontSize',14);ylabel('x','Interpreter','Latex','FontSize',14); title(sprintf('Solution to the advection equation using \n c=-(1+2sin(5t)-H(x-4))'),...
        'Interpreter','Latex'); colorbar();
hold on; yline(4,'r-','Linewidth',1.5)
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