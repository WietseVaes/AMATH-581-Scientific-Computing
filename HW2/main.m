close all; clear; clc;
%% Problem 1
L = 4; % max(|x|)
xp = linspace(-L,L,20*L+1); 
tol = 1e-6; %tol for break condition
beta_start = 0; % Initial beta guess
for modes = 1:5
    beta = beta_start; %starting beta for mode
    dbeta = 1/16; %starting change in beta
    for j = 1:10000
        A = sqrt(L.^2-beta); % Updating second initial condition y'(-L)
        y0 = [1; A];
        [x, y] = ode45( @(x, y) rhsfunc(x, y, beta), xp, y0); %Finding solution

        if abs(y(end, 2)+sqrt(L.^2-beta)*y(end,1))<tol %testing break condition : y'(L) approx ...
            break
        end

        if (-1)^(modes+1)*(y(end, 2)+sqrt(L.^2-beta)*y(end,1))>0 %Finding out how to change guess for beta
            beta = beta + dbeta;
        else
            beta = beta - dbeta/2;

            dbeta = dbeta/2;
        end
    end
    beta_start = beta + 0.1; %Changing starting beta for next mode.
    AA{modes} = y(:,1)/sqrt(trapz(xp,y(:,1).^2)); %Adding normalized Eigenvector to the group
    A6(modes) = beta; %set delivarable
end
%set delivarables
A1 = abs(AA{1});
A2 = abs(AA{2});
A3 = abs(AA{3});
A4 = abs(AA{4});
A5 = abs(AA{5});
%% Plotting the solution to the one-dimensional harmonic trapping potential using the second smallest eignvalue and respective eigenvector. 

[xx,tt] = meshgrid(x,linspace(0,5,1000)); % Making Mesh
t = linspace(0,5,1000);
Z = (AA{2}*cos(A6(2)*t/2))'; % Evaluating the psi_2 values
figure(); %creating a new figure
subplot(121);
surf(xx,tt,Z); shading interp %Plotting a 3D plot
xlabel('x','FontSize',17);ylabel('t','FontSize',17);zlabel('\psi_2(x,t)','FontSize',13); title('3D plot of \psi_2.','FontSize',13); %Labeling axis+title
subplot(122);
surf(xx,tt,Z); view(2)%Plotting a 3D plot
shading interp; colormap jet; h=colorbar; set(gcf,'renderer','zbuffer');%making colorbar
set(get(h,'title'),'string','\psi_2(x,t)','FontSize',12); xlabel('x','FontSize',17);ylabel('t','FontSize',17);title('Contour plot of \psi_2.','FontSize',13);%Labeling axis+colorbar+title
%}
%% Problem 2
L = 4;
x = linspace(-L,L,20*L+1);
x(1)=[];
x(end) = [];
Deltax = mean(diff(x));
n = length(x);
diago = -2-x.^2*Deltax^2;
sidediag = ones(1,n-1);
MatrixA = diag(diago) + diag(sidediag,1) + diag(sidediag,-1);
MatrixA(1,1) = MatrixA(1,1)+4/3;MatrixA(1,2) = MatrixA(1,2)-1/3;
MatrixA(n,n) = MatrixA(n,n)+4/3;MatrixA(n,n-1) = MatrixA(n,n-1)-1/3;
MatrixA = -MatrixA/(Deltax.^2);
%MatrixA = [zeros(size(MatrixA,1),1),MatrixA,zeros(size(MatrixA,1),1)]
%MatrixA = [1,zeros(1,size(MatrixA,2)-1);MatrixA;zeros(1,size(MatrixA,2)-1),1]
[V,E] = eig(MatrixA);
for index1 = 1:size(V,2)
    begpointvalue = (4*V(1,index1)-V(2,index1))/(3+sqrt(L^2-E(index1,index1))*2*Deltax);
    endpointvalue = (4*V(n,index1)-V(n-1,index1))/(3+sqrt(L^2-E(index1,index1))*2*Deltax);
    Eigenvector{index1} = [begpointvalue;V(:,index1);endpointvalue];
end
x = linspace(-L,L,20*L+1);

E(E==0) = []; Esorted = sort(E);
A12 = Esorted(1:5);
for index1 = 1:5
    indexsorted = find(E==Esorted(index1));
    CC{index1} = Eigenvector{indexsorted}/sqrt(trapz(x,Eigenvector{indexsorted}.^2));
end
%{
figure; hold on;
for index1 = 1:5
    plot(x,abs(CC{index1}))
end
%}
A7 = abs(CC{1});
A8 = abs(CC{2});
A9 = abs(CC{3});
A10 = abs(CC{4});
A11 = abs(CC{5});


clearvars -except A1 A2 A3 A4 A5 A6 A7 A8 A9 A10 A11 A12
%% Problem 3
n0 = 1;
L = 3;
xp = linspace(-L,L,20*L+1)'; % xpsan, don't define stepsize linspace(-L,L,20*L+1)
tol = 1e-5;
Gamma = [0.05,-0.05];
for gammaindex = 1:2
    gamma = Gamma(gammaindex);
    beta_start = 0.1;
    B = 10^(-3);
    for modes = 1:2
        beta = beta_start;
        dbeta = n0/16;
        for j = 1:10000
            A = sqrt(L.^2-beta)*B;
            y0 = [B; A];
            [x, y] = ode45( @(x, y) rhsfunc2(x, y, gamma, beta), xp, y0);

            if abs(y(end, 2)+sqrt(L.^2-beta)*y(end,1))<tol && abs(sqrt(trapz(xp,y(:,1).^2))-1) <= tol
                break
            else
               B = B/sqrt(trapz(xp,y(:,1).^2));
            end
            A = sqrt(L.^2-beta)*B;
            y0 = [B; A];
            [x, y] = ode45( @(x, y) rhsfunc2(x, y, gamma, beta), xp, y0);

            if abs(y(end, 2)+sqrt(L.^2-beta)*y(end,1))<tol && abs(sqrt(trapz(xp,y(:,1).^2))-1) <= tol
                break
            end
            if (-1)^(modes+1)*(y(end, 2)+sqrt(L.^2-beta)*y(end,1))>0
                beta = beta + dbeta; % Decrease beta, start at 100 and decrease
            else
                beta = beta - dbeta/2;

                dbeta = dbeta/2; % Redefine dbeta so we converge
                % faster / higher accuracy
            end
        end
        beta_start = beta + 0.1;
        BB{gammaindex}(:,modes) = y(:,1);
        Eigenvalues{gammaindex}(modes) = beta;
    end
end

%{
for index2 = 1:2
    figure; hold on
    for index1 = 1:2
        plot(x, BB{index2}(:,index1), 'linewidth', 2)
    end
end
%}
A13 = abs(BB{1}(:,1));
A14 = abs(BB{1}(:,2));
A15 = Eigenvalues{1};
A16 = abs(BB{2}(:,1));
A17 = abs(BB{2}(:,2));
A18 = Eigenvalues{2};


%% functions
function f = rhsfunc(t, y, beta)
f1 = y(2);
f2 = (t.^2 - beta).*y(1);
f = [f1; f2];
end
function f = rhsfunc2(t, y, gamma, beta)
f1 = y(2);
f2 = (gamma*y(1).^2+t.^2 - beta).*y(1);
f = [f1; f2];
end
%%
