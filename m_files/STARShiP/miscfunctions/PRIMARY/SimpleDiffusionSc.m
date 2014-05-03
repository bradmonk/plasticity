clc, clear all; close all;


[XL YL] = SimpleDiffusion();

XLC = reshape(XL,[],1);
YLC = reshape(YL,[],1);
XYLC = [XLC YLC];


XYLDat = [XYLC ones(numel(XLC),1) ones(numel(YLC),1)];

TIME = repmat((1:100)',10,1);
XYLDat(:,3) = TIME;

PART = reshape(repmat((1:10),100,1),[],1);
XYLDat(:,4) = PART;

XYLD = cat(1,[0 0 0 0],XYLDat);

xlswrite('DiffData',XYLD)







%---------------------------------------------------
% Export to JMP
%{.

% BASE DIFFUSION RATES EQUATIONS
Ns = 100;				% number of steps (loop time)
Np = 10;				% number of particles
D = 0.5;				% Diffusion Rate (D = L² / 2d*t)
dt = 1;					% time step (seconds)
dm = 2;                 % dimensions
k = sqrt(dm*D*dt);		% stdev of D's step size distribution
MSD = 2*dm*D;			% theoretical mean squared displacement
Tk = dm*k*(1:Ns);		% theoretical MSD over time
Lxy = zeros(2,Np);		% XY particle locations

% PERFORM BROWNIAN MOTION
for t = 1:Ns 
	Sxy = k * randn(2,Np);	% generates step sizes
	Lxy = Lxy+Sxy;			% adds step to location
	XL(t,:) = Lxy(1,:);		% saves x location
	YL(t,:) = Lxy(2,:);		% saves y location
end

% FORMAT MATRIX FOR JMP EXPORT
XLC = reshape(XL,[],1);
YLC = reshape(YL,[],1);
XYLC = [XLC YLC];
XYLDat = [XYLC ones(numel(XLC),1) ones(numel(YLC),1)];
TIME = repmat((1:100)',10,1);
XYLDat(:,3) = TIME;
PART = reshape(repmat((1:10),100,1),[],1);
XYLDat(:,4) = PART;
XYLD = cat(1,[0 0 0 0],XYLDat);
csvwrite('DiffData.csv',XYLD)
%}
%---------------------------------------------------



%---------------------------------------------------
% Brownian Motion from Standard Normal
%{
%---------------------------------------------------
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
% Einstein's Brownian Diffusion Gaussian Function
x = 0;		% origin
t = 250;	% time step
p0 = 100;	% particle density at origin
D = 0.1;	% mass diffusivity

xx=linspace(x-t*D,x+t*D,400);
pxt = (p0 ./ sqrt(4.*pi.*D.*t)) .* exp(-( (xx.^2) ./ (4.*D.*t)));
plot(xx,pxt);
%--------------------
%}
%---------------------------------------------------

