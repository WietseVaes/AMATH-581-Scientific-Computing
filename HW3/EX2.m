%% Problem 2
for m = [64,128]
    clearvars -except m TIME1 TIME2
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

    %Parameters
    nu = 0.001;
    f = @(x,y)  exp(-2*x.^2-(y.^2/20));
    L = 10;
    T = 4;
    Deltax = 20/m;
    Deltat = 0.5;

    xx = linspace(-L,L-Deltax,m)'; yy = xx;

    A = A/(Deltax^2);
    B = B/(2*Deltax);
    C = C/(2*Deltax);

    y = repmat(yy,m,1); %Creating x and y coordinates
    x = repmat(xx,1,m)';
    x=reshape(x,[m^2,1]);

    init = A\f(x,y); %initial condition
    [L,U,P] = lu(A); % LU decomposition

    %Solving without the LU decomposition
    tic
    [t, ppsi1] = ode45(@(t, x) vorticity1(t, x, A,B,C,nu), 0:Deltat:T, init);
    TIME1(m/64) = toc;

    %Solving with the LU decomposition
    tic
    [t, ppsi2] = ode45(@(t, x) vorticity2(t, x, A,B,C,nu,L,U,P), 0:Deltat:T, init);
    TIME2(m/64) = toc;
end
log(TIME1(2)/TIME1(1))/log(4)
log(TIME2(2)/TIME2(1))/log(4)


%% Functions
function ut = vorticity1(t, x, A,B,C,nu)
    b =  nu*A^2*x+(C*x).*(B*A*x)-(B*x).*(C*A*x);
    ut = A\b;
end
function ut = vorticity2(t, x, A,B,C,nu,L,U,P)
    b =  nu*A^2*x+(C*x).*(B*A*x)-(B*x).*(C*A*x);
    ut = U\(L\(P*b));
end