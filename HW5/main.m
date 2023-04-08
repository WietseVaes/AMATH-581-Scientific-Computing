clear all
%% problem 1
L = 20;
n = 64;
T = 25;
Deltat = .5;
tspan = 0:Deltat:T;
x = linspace(-L/2,L/2,n+1); x(end) = []; y=x;
[X, Y] = meshgrid(x,y);
A1 = X;
m = 3;
alpha = 0;
u0 = (tanh(sqrt(X.^2+Y.^2))-alpha).*cos(m*angle(X+1i*Y) - sqrt(X.^2+Y.^2));
v0 = (tanh(sqrt(X.^2+Y.^2))-alpha).*sin(m*angle(X+1i*Y) - sqrt(X.^2+Y.^2));

A2 = u0;

u0hat = fft2(u0);
A3 = real(u0hat);
v0hat = fft2(v0);
initcond = [reshape(u0hat,n^2,1); reshape(v0hat,n^2,1)];

A4 = imag(initcond);

kx = (2*pi/L)*[0:(n/2-1), (-n/2):-1]'; ky = kx;

[KX,KY] = meshgrid(kx,ky);
K = KX.^2+KY.^2;
Kvec = reshape(K,n^2,1);

Beta = 1;
D1 = .1;
D2 = .1;
[t,yvec] = ode45(@(t,yvec) FFT_reactdiff(t,yvec,K,Kvec,n,Beta,D1,D2),tspan,initcond);
yvec = yvec.';

for index1 = 1:length(t)
    Usol(:,:,index1) = real(ifft2(reshape(yvec(1:n^2,index1),n,n)));
    Vsol(:,:,index1) = real(ifft2(reshape(yvec((n^2+1):2*n^2,index1),n,n)));
end
%{
figure(3);
for index1 = 1:length(t)
    subplot(1,2,1)
    pcolor([Usol(:,:,index1)]); shading interp;
    subplot(1,2,2)
    pcolor([Vsol(:,:,index1)]); shading interp; drawnow;
end
%}
A5 = real(yvec);
A6 = imag(yvec);
A7 = real(yvec(1:n^2,5));
A8 = reshape(real(A7),n,n);
A9 = real(ifft2(reshape(yvec(1:n^2,5),n,n)));% figure(2);surf(X,Y,A9); shading interp

%% problem 2
clearvars -except A1 A2 A3 A4 A5 A6 A7 A8 A9
L = 20;
n = 30;
T = 25;
Deltat = .5;
tspan = 0:Deltat:T;

[D,x] = cheb(n);
D2 = D^2;
x = x*L/2;
D2 = 4/L^2*D2;

D2 = D2(2:n,2:n);
x2 = x(2:n); y=x2;
[X, Y] = meshgrid(x2,y);
A11 = Y;
m = 2;
alpha = 1;
u0 = (tanh(sqrt(X.^2+Y.^2))-alpha).*cos(m*angle(X+1i*Y) - sqrt(X.^2+Y.^2));
v0 = (tanh(sqrt(X.^2+Y.^2))-alpha).*sin(m*angle(X+1i*Y) - sqrt(X.^2+Y.^2));

A12 = v0; 

initcond = [reshape(u0,(n-1)^2,1); reshape(v0,(n-1)^2,1)];

A13 = initcond;

I = eye(n-1);

Lap = kron(D2,I)+ kron(I,D2);

A10 = Lap;

Beta = 1;
D1 = .1;
D2 = .1;
[t, ysolvec] = ode45(@(t,yvec) FFT_reactdiff_cheb(t,yvec,Lap,n-1,Beta,D1,D2), tspan, initcond);
A14 = ysolvec;

A15 = ysolvec(5, (n-1)^2+1:2*(n-1)^2)';

for index1 = 1:length(t)
    Usol(:,:,index1) = reshape(ysolvec(index1,1:(n-1)^2)',n-1,n-1);
    Vsol(:,:,index1) = reshape(ysolvec(index1,((n-1)^2+1):2*(n-1)^2)',n-1,n-1);

    U(:,:,index1) = [zeros(n+1,1),[zeros(1,n-1);Usol(:,:,index1);zeros(1,n-1)],zeros(n+1,1)];
    V(:,:,index1) = [zeros(n+1,1),[zeros(1,n-1);Vsol(:,:,index1);zeros(1,n-1)],zeros(n+1,1)];
end
A16 = V(:,:,5);
%{
figure(4);
for index1 = 1:length(t)
    subplot(1,2,1)
    pcolor(U(:,:,index1)); shading interp;
    subplot(1,2,2)
    pcolor(V(:,:,index1)); shading interp; drawnow;
end
%}
%% functions
function rhs = FFT_reactdiff(t,yvec,K,Kvec,n,beta,D1,D2)
    uffvec = yvec(1:n^2);
    vffvec = yvec((n^2+1):2*n^2);
    uff = reshape(uffvec,n,n);
    vff = reshape(vffvec,n,n);
    u = real(ifft2(uff));
    v = real(ifft2(vff));
    Asqrd = abs(u).^2+abs(v).^2;
    lambda = 1-Asqrd;
    omega = -beta*Asqrd;
    rhs(1:n^2,1) = reshape(fft2(lambda.*u - omega.*v)- D1*K.*uff,n^2,1) ;
    rhs((n^2+1):2*n^2,1) = reshape(fft2(omega.*u-lambda.*v)- D2*K.*vff ,n^2,1) ;
end

function rhs = FFT_reactdiff_cheb(t,yvec,Lap,n,beta,D1,D2)
    uvec = yvec(1:n^2);
    vvec = yvec((n^2+1):2*n^2);
    u = reshape(uvec,n,n);
    v = reshape(vvec,n,n);
    
    Asqrd = abs(u).^2+abs(v).^2;
    lambda = 1-Asqrd;
    omega = -beta*Asqrd;
    rhs(1:n^2,1) = reshape(lambda.*u - omega.*v,n^2,1)+ D1*Lap*uvec;
    rhs((n^2+1):2*n^2,1) = reshape(-lambda.*v + omega.*u,n^2,1)+ D2*Lap*vvec ;
end
