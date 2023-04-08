
B0 = 1;
gamma = 1/10;
W = 3;
t = linspace(-10,10,1000);
y0 = [0,B0];

[t, y] = ode45( @(t, y) fuckme(t,y), t, y0);
figure(1);
plot(t,y(:,1))
figure(2);
plot(t,y(:,2))

function dfdt=fuckme(t,y)
B0 = 2;
gamma = 1/10;
W = 3;
    f1 = y(2).*(y(2).^2+y(1).^2-W)-B0*(B0^2-W);
    f2 = -y(1).*(y(2).^2+y(1).^2-gamma*B0^2-W);
    dfdt = [f1; f2];
end

