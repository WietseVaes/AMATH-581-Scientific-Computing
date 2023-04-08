clear all; close all; clc;

%% Problem 1
func = @(t) pi*exp(3*(cos(t)-1))/sqrt(2);
f = @(t,x) -3*x*sin(t);
DeltaT = 2.^(-(2:8));
for i = 1:length(DeltaT)
    clear y
    y(1,1:3) = pi/sqrt(2);
    dt = DeltaT(i);
    t = 0:dt:5;
    for j = 1:(length(t)-1)
        y(j+1,1) = y(j,1) + dt*f(t(j),y(j,1));
        y(j+1,2) = y(j,2) + dt/2*(f(t(j),y(j,2))+f(t(j+1),y(j,2)+dt*f(t(j),y(j,2))));
        if j ~= 1
            yp = y(j,3) + dt/2*(3*f(t(j),y(j,3))-f(t(j-1),y(j-1,3)));
            y(j+1,3) = y(j,3) + dt/2*(f(t(j+1),yp)+f(t(j),y(j,3)));
        else
            y(j+1,3) = y(3)+ dt*f(t(j)+dt/2,y(j,3)+dt/2*f(t(j),y(j,3)));
        end
    end
    for j = 1:3
        E(j,i) = abs(func(5)-y(end,j));
    end
end
for j = 1:3
    p(j,:) = polyfit(log(DeltaT),log(E(j,:)),1);
end
%{
loglog(DeltaT,E(1,:),"*-")
hold on
plot(DeltaT,(DeltaT-DeltaT(end))+E(1,end))
figure();
loglog(DeltaT,E(2,:),"*-")
hold on
loglog(DeltaT,(DeltaT-DeltaT(end)).^2+E(2,end))
figure();
loglog(DeltaT,E(3,:),"*-")
hold on
loglog(DeltaT,(DeltaT-DeltaT(end)).^2+E(3,end))
for j = 1:3
    yy(j,:) = polyval(p(j,:),log(DeltaT));
end
%loglog(DeltaT,yy')
legend("Forward-Euler", "Heun","Adams")
hold off
%}
A1 = y(:,1);
A2 = E(1,:);
A3 = p(1,1);
A4 = y(:,2);
A5 = E(2,:);
A6 = p(2,1);
A7 = y(:,3);
A8 = E(3,:);
A9 = p(3,1);



%% Problem 2
%a
clear y t p
Epsilon = [0.1,1,20];
y0 = [sqrt(3);1];
tspan = 0:0.5:32;
for i = 1:3
    epsilon = Epsilon(i);
    f = @(t,y) [y(2); -(epsilon*(y(1)^2-1)*y(2)+y(1))];
    [t,ysol] = ode45(f, tspan, y0);
    A10(:,i) = ysol(:,1);
end

%plot(t,ysol(:,1),'-o')

%b
f = @(t,y) [y(2); -((y(1)^2-1)*y(2)+y(1))];
y0 = [2;pi^2];
Tols = 10.^-(4:10);
for i = 1:7
    clear options tol T1 T2 T3
    tol = Tols(i);
    options = odeset('AbsTol', tol, 'reltol', tol);
    [T1, Y1] = ode45(f, [0,32], y0, options);
    DeltaTavg(i,1) = mean(diff(T1));

    [T2, Y2] = ode23(f, [0,32], y0, options);
    DeltaTavg(i,2) = mean(diff(T2));

    [T3, Y3] = ode15s(f, [0,32], y0, options);
    DeltaTavg(i,3) = mean(diff(T3));
end
p = polyfit(log(DeltaTavg(:,1)),log(Tols),1);
A11 = p(1);
p = polyfit(log(DeltaTavg(:,2)),log(Tols),1);
A12 = p(1);
p = polyfit(log(DeltaTavg(:,3)),log(Tols),1);
A13 = p(1);

%% Problem 3
clear times T Y a1 a2 b c I f y0
a1 = 0.05;
a2 = 0.25;
b = 0.1;
c = 0.1;
I = 0.1;
DD = [0,0;0,0.2;-0.1,0.2;-0.3,0.2;-0.5,0.2;0.1,0.2];
y0 = [0.1;0.1;0;0];
for i = 1:size(DD,1)
    d = DD(i,:);
    f = @(t,y) [-y(1).^3+(1+a1).*y(1).^2-a1*y(1)-y(3)+I+d(1).*y(2); -y(2).^3+(1+a2)*y(2).^2-a2*y(2)-y(4)+I+d(2)*y(1);b*y(1)-c*y(3);b*y(2)-c*y(4)];
    [T, Y(:,:,i)] = ode15s(f, 0:0.5:100, y0);
end
A14 = Y(:,:,1);
A15 = Y(:,:,2);
A16 = Y(:,:,3);
A17 = Y(:,:,4);
A18 = Y(:,:,5);
A19 = Y(:,:,6);
%
plot(T,A14');
legend('v1','v2','w1','w2')
xlabel('t')
title('FH model approximation with interaction parameter (0,0)')
figure();
subplot(2,2,1)
plot(T,A15');
legend('v1','v2','w1','w2')
xlabel('t')
title('Interaction parameter (0,0.2)')
subplot(2,2,2)
plot(T,A16');
legend('v1','v2','w1','w2')
xlabel('t')
title('Interaction parameter (-0.1,0.2)')
subplot(2,2,3)
plot(T,A17');
legend('v1','v2','w1','w2')
xlabel('t')
title('Interaction parameter (-0.3,0.2)')
subplot(2,2,4)
plot(T,A18');
legend('v1','v2','w1','w2')
title('Interaction parameter (-0.5,0.2)')
xlabel('t')
%}
figure();
plot(T,A19');
legend('v1','v2','w1','w2')
title('Interaction parameter (0.1,0.2)')
xlabel('t')

