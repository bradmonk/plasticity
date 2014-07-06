
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
axis equal;
%axis([-5 5 0 1.5])
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
% 3D Gaussian Distribution
A = 5;
x0 = 0; y0 = 0;
sigma_x = 1;
sigma_y = 3;

a = cos(theta)^2/2/sigma_x^2 + sin(theta)^2/2/sigma_y^2;
b = -sin(2*theta)/4/sigma_x^2 + sin(2*theta)/4/sigma_y^2 ;
c = sin(theta)^2/2/sigma_x^2 + cos(theta)^2/2/sigma_y^2;

[X, Y] = meshgrid(-5:.1:5, -5:.1:5);
Z = A*exp( - (a*(X-x0).^2 + 2*b*(X-x0).*(Y-y0) + c*(Y-y0).^2)) ;

surf(X,Y,Z);
shading interp; view(-36,10);axis equal;

 
for theta = 0:pi/100:pi
    a = cos(theta)^2/2/sigma_x^2 + sin(theta)^2/2/sigma_y^2;
    b = -sin(2*theta)/4/sigma_x^2 + sin(2*theta)/4/sigma_y^2 ;
    c = sin(theta)^2/2/sigma_x^2 + cos(theta)^2/2/sigma_y^2;
 
    [X, Y] = meshgrid(-5:.1:5, -5:.1:5);
    Z = A*exp( - (a*(X-x0).^2 + 2*b*(X-x0).*(Y-y0) + c*(Y-y0).^2)) ;
    surf(X,Y,Z);shading interp;view(-36,10);
	axis equal;
	%axis vis3d;
	drawnow
end
%--------------------


%--------------------
A  = 1;		% hight of peak
x0 = 0;		% x-axis peak location 
y0 = 0;		% y-axis peak location 
sd = .5;	% sigma (stdev of slope)
tta = 0;	% theta (axis revolution)
res= .2;	% grid resolution
spr= 1;		% grid spread	

a = cos(tta)^2/2/sd^2 + sin(tta)^2/2/sd^2;
b = -sin(2*tta)/4/sd^2 + sin(2*tta)/4/sd^2 ;
c = sin(tta)^2/2/sd^2 + cos(tta)^2/2/sd^2;

[X, Y] = meshgrid((-spr):(res):(spr), (-spr):(res):(spr));
Z = A*exp( - (a*(X-x0).^2 + 2*b*(X-x0).*(Y-y0) + c*(Y-y0).^2)) ;

surf(X,Y,Z);
	view(-45,30); 
	axis([-spr spr -spr spr 0 A])
	%axis equal; 
	%axis vis3d; 
	shading interp;
	xlabel('x-axis');ylabel('y-axis');zlabel('z-axis')
%--------------------
	
	

%--------------------
A  = 1;		% hight of peak
x0=0;y0=0;	% x&y-axis peak locations
sd = .2;	% sigma (stdev of slope)
a = .5/sd^2;
c = .5/sd^2;

% Gres = .1;				% grid resolution
% Gspr = 1;					% grid spread	
% Gnum = 2*Gspr/Gres +1;	% grid N (odd num)
Gres = .1;
Gnum = 11;				
Gspr = ((Gnum-1)*Gres)/2;

[X, Y] = meshgrid((-Gspr):(Gres):(Gspr), (-Gspr):(Gres):(Gspr));
Z = A*exp( - (a*(X-x0).^2 + c*(Y-y0).^2)) ;

%---
surf(X,Y,Z);
view(-45,30); 
axis([-Gspr Gspr -Gspr Gspr 0 A])
shading interp; %axis equal; %axis vis3d; 
xlabel('x-axis');ylabel('y-axis');zlabel('z-axis')
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
% Change Factory Presents (unrelated)
a = get(0,'Factory');
set(0,'DefaultAxesFontName','Century Gothic')
set(0,'DefaultAxesFontSize',12)
factoryAxesPosition: [0.1300 0.1100 0.7750 0.8150]
factoryAxesTickLength: [0.0100 0.0250]
%--------------------















