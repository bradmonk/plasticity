clc, close all; clear all; scsz = get(0,'ScreenSize');

%---------------------------------------------------------------
% 3D Gaussian Distribution
fh5 = figure(5); set(fh5,'OuterPosition',(scsz./[2e-3 2e-3 2 2]))
%---------------------------------------------------------------
% A = 2;		res=.1;	x0=(-res/2); y0=(-res/2);	sx = .2; sy = .2;	
A = 2;		res=.2;	x0=0; y0=0;	sx = .2; sy = .2;	

t = 0;
a = cos(t)^2/2/sx^2 + sin(t)^2/2/sy^2;
b = -sin(2*t)/4/sx^2 + sin(2*t)/4/sy^2 ;
c = sin(t)^2/2/sx^2 + cos(t)^2/2/sy^2;

[X, Y] = meshgrid((-sx*3):(res):(sx*3), (-sy*3):(res):(sy*3));
Z = A*exp( - (a*(X-x0).^2 + 2*b*(X-x0).*(Y-y0) + c*(Y-y0).^2)) ;

figure(fh5)
subplot('Position',[.05 .05 .40 .90]); ph5 = surf(X,Y,Z);
view(0,90); axis equal; drawnow;
xlabel('x-axis');ylabel('y-axis');zlabel('z-axis')
subplot('Position',[.55 .05 .40 .90]); ph6 = surf(X,Y,Z);
view(-13,22); % shading interp; 
axis equal;drawnow;
xlabel('x-axis');ylabel('y-axis');zlabel('z-axis')
%----------------------
for t = 0:pi/100:pi
	a = cos(t)^2/2/sx^2 + sin(t)^2/2/sy^2;
	b = -sin(2*t)/4/sx^2 + sin(2*t)/4/sy^2 ;
	c = sin(t)^2/2/sx^2 + cos(t)^2/2/sy^2;

	[X, Y] = meshgrid((-sx*3):(res):(sx*3), (-sy*3):(res):(sy*3));
	Z = A*exp( - (a*(X-x0).^2 + 2*b*(X-x0).*(Y-y0) + c*(Y-y0).^2)) ;

	set(ph5,'XData',X,'YData',Y,'ZData',Z);
	drawnow;
	set(ph6,'XData',X,'YData',Y,'ZData',Z);
	drawnow;
end
%---------------------------------------------------------------


zmask=	[0.0001    0.0036    0.0111    0.0036    0.0001;
		 0.0036    0.1054    0.3247    0.1054    0.0036;
		 0.0111    0.3247    1.0000    0.3247    0.0111;
		 0.0036    0.1054    0.3247    0.1054    0.0036;
		 0.0001    0.0036    0.0111    0.0036    0.0001];

vals = SAPMASKS(LBR,TIME,SIZE,MODS,DOES,REVA,GLU,GT,GTab)




%--------------------
% 3D Gaussian Distribution
A = 5;
x0 = 0; y0 = 0;
sx = 1;
sy = 3;
 
for t = 0:pi/100:pi
    a = cos(t)^2/2/sx^2 + sin(t)^2/2/sy^2;
    b = -sin(2*t)/4/sx^2 + sin(2*t)/4/sy^2 ;
    c = sin(t)^2/2/sx^2 + cos(t)^2/2/sy^2;
 
    [X, Y] = meshgrid(-5:.1:5, -5:.1:5);
    Z = A*exp( - (a*(X-x0).^2 + 2*b*(X-x0).*(Y-y0) + c*(Y-y0).^2)) ;
    surf(X,Y,Z);shading interp;view(-36,10);axis equal;drawnow
end
%--------------------







%--------------------
% Standard Normal Distribution
mu=0;
sd=1;
x=linspace(mu-4*sd,mu+4*sd,200);
pdfx=1/sqrt(2*pi)/sd*exp(-(x-mu).^2/(2*sd^2));
plot(x,pdfx);
%--------------------



%--------------------
% Universal Gaussian Function
a = 4;	% hight of peak from baseline (Yaxis)
b = 0;	% mean (peak center on Xaxis)
c = 1;	% stdev (curve spread)
d = 0;	% baseline (Yaxis)
x = -4:.1:4;

fx = a .* exp(-((x-b).^2 ./ (2*c).^2)) + d;

plot(x,fx);
axis([-5 5 0 1.5])
hold on
%--------------------
% Universal Gaussian Function
a = 10;	% hight of peak from baseline (Yaxis)
b = 0;	% mean (peak center on Xaxis)
c = 1;	% stdev (curve spread)
d = 0;	% baseline (Yaxis)
x = (b-5*c):.1:(b+5*c);

fx = a .* exp(-((x-b).^2 ./ (2*c).^2)) + d;

plot(x,fx,'r');
axis([(b-5*c) (b+5*c) d ((a+d)*1.2)])
%--------------------



%--------------------
% Gaussian Function Dual
a = 10;	% hight of peak from baseline (Yaxis)
b = 20;	% mean (peak center on Xaxis)
c = 2;	% stdev (curve spread)
d = 0;	% baseline (Yaxis)
x = (b-5*c):.1:(b+5*c);

fx = a .* exp(-((x-b).^2 ./ (2*c).^2)) + d;

figure(1)
plot(x,fx,'r.');
axis([(0) (80) d ((a+d)*1.5)])
hold on

% Gaussian Function
a = 10;	% hight of peak from baseline (Yaxis)
b = 40;	% mean (peak center on Xaxis)
c = 4;	% stdev (curve spread)
d = 0;	% baseline (Yaxis)
x = (b-5*c):.1:(b+5*c);

fx = a .* exp(-((x-b).^2 ./ (2*c).^2)) + d;

figure(1)
plot(x,fx,'b.');
axis([(0) (80) d ((a+d)*1.5)])
xlabel('seconds'); ylabel('responses')

%--------------------



%--------------------
% Einstein's Brownian Diffusion Gaussian Function
x = 0;		% origin
t = 250;	% time step
p0 = 100;	% particle density
D = 0.1;	% mass diffusivity

xx=linspace(x-t*D,x+t*D,400);
pxt = (p0 ./ sqrt(4.*pi.*D.*t)) .* exp(-( (xx.^2) ./ (4.*D.*t)));
plot(xx,pxt);
%--------------------




%--------------------
% Dirac Delta Function
for a = 1:10
x=linspace(-20,20,100);
dax = (1 ./ (a .* sqrt(pi)) ) .* exp(-x.^2 ./ a.^2);
figure(10)
plot(x,dax);
hold on
pause(1)
end
%--------------------


%---------------------------------------------------
% Brownian Motion from Standard Normal
N = 1000;
displacement = randn(1,N);
plot(displacement);
particle = struct();
particle.x = cumsum( randn(N, 1) );
particle.y = cumsum( randn(N, 1) );
plot(particle.x, particle.y);
ylabel('Y Position');
xlabel('X Position');
title('position versus time in 2D');
%--------------------
dsquared = particle.x .^ 2 + particle.y .^ 2;
plot(dsquared);
%--------------------
D = 0.1;				% Diffusion coefficient um^2/s
dim = 2;	         % two dimensional simulation
tau = .1;               % time interval in seconds
time = tau * 1:N;       % create a time vector for plotting
k = sqrt(D * dim * tau);
dx = k * randn(N,1);
dy = k * randn(N,1);
x = cumsum(dx);
y = cumsum(dy);
dSquaredDisplacement = (dx .^ 2) + (dy .^ 2);
squaredDisplacement = ( x .^ 2) + ( y .^ 2);
plot(x,y);
title('Particle Track of a Single Simulated Particle');
%--------------------
clf;
hold on;
plot(time, (0:1:(N-1)) * 2*k^2 , 'k', 'LineWidth', 3);      % plot theoretical line
plot(time, squaredDisplacement);
hold off;
xlabel('Time');
ylabel('Displacement Squared');
title('Displacement Squared versus Time for 1 Particle in 2 Dimensions');
%--------------------
simulatedD = mean( dSquaredDisplacement ) / ( 2 * dim * tau );
standardError = std( dSquaredDisplacement ) / ( 2 * dim * tau * sqrt(N) );
actualError = D - simulatedD;
disp([simulatedD standardError actualError])
%---------------------------------------------------


%--------------------
% 3D Gaussian Distribution
A = 5;
x0 = 0; y0 = 0;
sx = 1;
sy = 3;
 
for t = 0:pi/100:pi
    a = cos(t)^2/2/sx^2 + sin(t)^2/2/sy^2;
    b = -sin(2*t)/4/sx^2 + sin(2*t)/4/sy^2 ;
    c = sin(t)^2/2/sx^2 + cos(t)^2/2/sy^2;
 
    [X, Y] = meshgrid(-5:.1:5, -5:.1:5);
    Z = A*exp( - (a*(X-x0).^2 + 2*b*(X-x0).*(Y-y0) + c*(Y-y0).^2)) ;
    surf(X,Y,Z);shading interp;view(-36,10);axis equal;drawnow
end
%--------------------




%--------------------
% Change Factory Presents (unrelated)
a = get(0,'Factory');
set(0,'DefaultAxesFontName','Century Gothic')
set(0,'DefaultAxesFontSize',12)
factoryAxesPosition: [0.1300 0.1100 0.7750 0.8150]
factoryAxesTickLength: [0.0100 0.0250]
%--------------------

















