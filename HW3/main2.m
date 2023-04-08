clear all
%% Problem 2
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
Deltat = 0.01;

xx = linspace(-L,L-Deltax,64)'; yy = xx;

A = A/(Deltax^2);
B = B/(2*Deltax);
C = C/(2*Deltax);

y = repmat(yy,64,1);
x = repmat(xx,1,64)';
x=reshape(x,[64^2,1]);

init = A\f(x,y);
[L,U,P] = lu(A);

tic
[t, ppsi] = ode45(@(t, x) vorticity2(t, x, A,B,C,nu,L,U,P), 0:Deltat:T, init);
toc

for index1 = 1:length(t)
    Omega(index1,:) = A*ppsi(index1,:)';

end
for index1 = 1:length(t)
    omega(:,:,index1) = reshape(Omega(index1,:),[64,64]);
    psi(:,:,index1) = reshape(ppsi(index1,:)',[64,64]);
end
omega(omega(:,:,:)<=0)=0;
omega(:,65,:) = zeros(64,1,length(t));omega(65,:,:) = zeros(1,65,length(t));
[XX,YY] = meshgrid([xx;10],[yy;10]);
%% animation

filename = 'heat_animation.gif';

figure(1);
[M,p] = contourf(XX,YY,omega(:,:,1),100,'edgecolor','none');colormap jet;
 h=colorbar; h.Limits = [0,.99]; xlabel('x','Interpreter','Latex','FontSize',14); ylabel('y','Interpreter','Latex','FontSize',14)
hold on
title(sprintf('Time evolution of the vorticity \x03C9(x, y, t) \n Time t = %0.2f', 0),...
        'Interpreter','Latex','Fontsize',10);
set(get(gca,'title'),'Position',[0 9.95 1.00011])

for k = 1:length(t)
    p.ZData = omega(:,:,k);
    title(sprintf('Time evolution of the vorticity \x03C9(x, y, t) \n Time t = %0.2f', t(k)),...
        'Interpreter','Latex');
    
    frame = getframe(gcf);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if k ==1
        imwrite(imind,cm,filename,'gif','LoopCount',inf,...
            'DelayTime',0.001)
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append',...
            'DelayTime',0.001)
    end
    
end

%% Functions
function ut = vorticity2(t, x, A,B,C,nu,L,U,P)
    b =  nu*A^2*x+(C*x).*(B*A*x)-(B*x).*(C*A*x);
    ut = U\(L\(P*b));
end