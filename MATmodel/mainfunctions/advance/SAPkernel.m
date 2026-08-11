function [hkMask] = SAPkernel(MSK)
% Creates a kernel mask using a 3D Gaussian
% that can be used for SAP clustering and
% other matrix conv() functions


GNpk = MSK{1};	% hight of peak
GNx0 = MSK{2};	% x-axis peak locations
GNy0 = MSK{3};	% y-axis peak locations
GNsd = MSK{4};	% sigma (stdev of slope)

GNnum = MSK{5};
GNres = MSK{6};
GNspr = ((GNnum-1)*GNres)/2;

a = .5/GNsd^2;
c = .5/GNsd^2;

[X, Y] = meshgrid((-GNspr):(GNres):(GNspr), (-GNspr):(GNres):(GNspr));
Z = GNpk*exp( - (a*(X-GNx0).^2 + c*(Y-GNy0).^2)) ;


hkMask=Z;

% FIGURE: 3D Gaussian Distribution
fh5 = figure(5);
set(fh5,'OuterPosition',[200 200 600 500]);
%----------------%
    figure(fh5)
subplot('Position',[.05 .55 .30 .40]); 
	ph5 = imagesc(hkMask); 
	axis equal;
subplot('Position',[.04 .08 .32 .42]); 
	ph7 = surf(X,Y,Z);
	axis equal; shading interp; view(90,90); 
subplot('Position',[.45 .05 .50 .90]); 
	ph7 = surf(X,Y,Z);
	axis vis3d; shading interp;
	view(-45,30); 
	xlabel('x-axis');ylabel('y-axis');zlabel('z-axis')


end
