clear; close all; clc

%% Introducing the FFT                      

L = 20; % define the computational domain [-L/2,L/2]
n = 128; % define the number of Fourier modes 2^n

x2 = linspace(-L/2,L/2,n+1); % define the domain discretization
x = x2(1:n);% consider only the first n points: periodicity
%(first and last points are the same)
u = exp(-x.*x); %function to take the fft of

ut = fft(u); %FFT the function

utshift = fftshift(ut); %shift FFT

%mention k space shift

figure(1), plot(x,u) %plot initial gaussian
figure(2), plot(abs(ut)) % plot unshifted transform
figure(3), plot(abs(utshift)) % plot shifted transform

%% Using FFT to take derivatives

L = 20; % define the computational domain [-L/2,L/2]
n = 128; % define the number of Fourier modes 2^n

x2 = linspace(-L/2,L/2,n+1); % define the domain discretization
x = x2(1:n);% consider only the first n points: periodicity

u = sech(x); %function to take the derivative of 
ut = fft(u); % FFT of the function

k = (2*pi/L)*[0:(n/2-1) (-n/2):-1]; % k rescaled to 2pi domain
% we don't have to do fft here, we just need to remember that everything is
% rotated
% to get the derivative multiply ut by i*k, repeat for more 
ut1 = i*k.*ut; % first derivative
ut2 = -k.*k.*ut; % second derivative
ut3 = -i*k.*k.*k.*ut; % third derivative


u1 = ifft(ut1); %inverse transform
u2 = ifft(ut2); %inverse transform
u3 = ifft(ut3); %inverse transform

u1exact = -sech(x).*tanh(x); %analytic first derivative

figure()
plot(x,u,'r')
hold on
plot(x,u1,'g')
plot(x,u1exact,'go')
plot(x,u2,'b')
plot(x,u3,'c')
legend('funcition','d1 approx','d1 true','d2 approx','d3 approx')

% if not smooth, increase points but remember powers of 2, so 256

%% Solve the Heat equation using FFT
% u_t = 5 u_{xx}, u(x,0) = f(x),
% |u| -> 0 as |x| -> infinity
% Use periodic boundary conditions and big L

f = @(x) exp(-x.^2);
L = 20;
n = 128;
T = 2;
t = 0:0.01:T;

x2= linspace(-L/2,L/2,n+1);
x = x2(1:n);

y0 = f(x);

y0hat = fft(y0);

k = (2*pi/L)*[0:(n/2-1) (-n/2):-1].';

FFT_heat = @(t,u) -5*k.*k.*u;

tic
[t,yhat] = ode45(FFT_heat,t,y0hat);
toc

y = real(ifft(yhat,n,2));

figure()
[X,T] = meshgrid(x,t);
surf(X,T,y,'DisplayName','Temperature Surface')
shading interp
xlabel('x,space')
ylabel('t,time')
zlabel('Temperature')

hold on

plot3(x,0*x,f(x),'r','linewidth',3,'DisplayName','Initial Heat Distribution')

legend()
title('Temperature modeled by the heat equation')

%% contour plot

figure()
contourf(X,T,y,'edgecolor','none')
colorbar()

%%

figure()
s = pcolor(X,T,y);
set(s,'EdgeColor','none')
s.FaceColor = 'interp';
colorbar()

%% animation

filename = 'heat_animation.gif';

figure
plot(x,y(1,:),'Color','none')
xlim([x(1),x(end)])
ylim([0,1])
hold on
p=plot(x,y(1,:),'b','linewidth',2)

for k = 1:length(t)
    p.XData = x;
    p.YData = y(k,:);
    title(sprintf('Heat solution\n Time: %0.2f sec', t(k)),...
        'Interpreter','Latex');
    
    pause(0.01)
    
    frame = getframe(gcf);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if k ==1
        imwrite(imind,cm,filename,'gif','LoopCount',inf,...
            'DelayTime',0.1)
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append',...
            'DelayTime',0.1)
    end
    
end

%%

figure();

h = figure;
plot(x,y(1,:),'linewidth',2)
axis tight manual
xlim([x(1),x(end)])
ylim([0,1])
xlabel('x')
ylabel('Magnitude of heat')
hold on
set(gca,'nextplot','replacechildren')

v = VideoWriter('heat_animation.mp4','MPEG-4')
open(v)


for j = 1:length(t)
    plot(x,y(j,:),'linewidth',2);
    xlim([x(1),x(end)])
    ylim([0,1])
    title(sprintf('Heat solution\n Time: %0.2f sec', t(j)),...
        'Interpreter','Latex');
    
    xlabel('x')
    ylabel('Magnitude of heat')
    
    frame = getframe(gcf);
    writeVideo(v,frame);
end
close(v);