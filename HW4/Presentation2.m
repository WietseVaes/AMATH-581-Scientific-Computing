clear all
%% Parameters + Getting Lambda
L = 10;
T = 2;
alpha = 2;
nt = 501;

x = linspace(-L,L,128+1); x(end)=[];
Deltax = mean(diff(x));
t = linspace(0,T,nt); Deltat = mean(diff(t));

lambda = Deltat*alpha/(Deltax^2);


%% calculating approximations
%using bicgstab
[u3,x,t] = calcubicgstab(512,lambda);

%Making u3 plotable.
u3(size(u3,1)+1,:) = u3(1,:);
x(length(x)+1) = 10;
nt = length(t);
%% plotting
[TT,XX] = meshgrid(t,x);
clf;
figure(1);
TTT = [1,nt,1000,2000,4000]; %Timestamp when solution will be explicitly shown
angle = [45,45;-45,45]; %angles of plotting
for index2 = 1:2
    subplot(1,2,index2);hold on
    surf(TT,XX,u3); shading interp;colormap turbo; %Plotting full solution
    for index1 = 1:length(TTT)
        plot3(ones(length(x),1)*t(TTT(index1)),x,u3(:,TTT(index1)),'Color','r','LineWidth',2); %plotting solution on specific timestamps
    end
    view(angle(index2,:));
    xlabel('t','Interpreter','Latex','FontSize',16);ylabel('x','Interpreter','Latex','FontSize',16);zlabel('u(x,t)','Interpreter','Latex','FontSize',16);
    hold off
end
sgtitle('Heat solution u(x,t) over time t.','Interpreter','Latex','FontSize',19)
%% Function
function [u, x, t] = calcubicgstab(nx,lambda)
%parameters
L = 10;
T = 2;
alpha = 2;

%grid definition
x = linspace(-L,L,nx+1); x(end)=[];
Deltax = mean(diff(x));
Deltat = lambda/alpha*Deltax^2;

t = 0:Deltat:T; nt = length(t);

%making matrix B and C
e1 = ones(nx,1); 
B = spdiags([-lambda*e1/2,e1,-lambda*e1/2],-1:1, nx, nx); C = spdiags([lambda*e1/2, e1, lambda*e1/2],-1:1, nx, nx);
B = B + lambda*speye(nx,nx); C = C - lambda*speye(nx,nx);
B(1,end) = -lambda/2;B(end,1) = -lambda/2;
C(1,end) = lambda/2;C(end,1) = lambda/2;

%initial condition
f = @(x) 10*cos(2*pi*x/L)+30*cos(8*pi*x/L);
u(:,1) = f(x)';
for index1 = 2:nt
    [u(:,index1),~] = bicgstab(B,C*u(:,index1-1)); %computing solution
end

end