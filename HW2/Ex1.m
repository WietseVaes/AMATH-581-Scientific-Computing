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
%% Plotting the solution to the one-dimensional harmonic trapping potential using the second smallest eigenvalue and respective eigenvector. 

[xx,tt] = meshgrid(x,linspace(0,5,1000)); % Making Mesh
t = linspace(0,5,1000);
Z = (AA{2}*cos(A6(2)*t/2))'; % Evaluating the psi_2 values
figure(); %creating a new figure
subplot(121);
surf(xx,tt,Z); shading interp %Plotting a 3D plot
xlabel('x','FontSize',17);ylabel('t','FontSize',17);zlabel('\psi_2(x,t)','FontSize',15); title('3D plot of \psi_2.','FontSize',15); %Labeling axis+title
subplot(122);
surf(xx,tt,Z); view(2)%Plotting a 3D plot
shading interp; colormap jet; h=colorbar; set(gcf,'renderer','zbuffer');%making colorbar
set(get(h,'title'),'string','\psi_2(x,t)','FontSize',12); xlabel('x','FontSize',17);ylabel('t','FontSize',17);title('Contour plot of \psi_2.','FontSize',15);%Labeling axis+colorbar+title

%% functions
function f = rhsfunc(t, y, beta)
f1 = y(2);
f2 = (t.^2 - beta).*y(1);
f = [f1; f2];
end