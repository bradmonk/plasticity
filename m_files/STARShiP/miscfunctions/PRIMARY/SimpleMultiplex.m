function [varargout] = SimpleMultiplex(varargin)
format compact;format short;close all;clc; clear all;
scsz = get(0,'ScreenSize');


% USER ENTERED VALUES
NSteps = 5000;			% number of steps (loop time)
Ndots = 100;			% number of particles
% Scale = 1/10;			% scale of model (life:model)
Scale = 1;				% scale of model (life:model)
TimeStep = 1000;		% time step (ms)
Delay = 0;				% slow animation rate by X%
POLYSz = [100 100];		% size of polygon enclosure (XY in µm)
DiffRateA = 4.0;		% diffusion rate coefficient A
DiffRateB = 0.1;		% diffusion rate coefficient B


% BASE DIFFUSION RATES EQUATIONS
Sc = Scale;					% scale of model (life:model)
t = TimeStep/1000;			% time step (ms)
dm = 2;                     % dimensions
Da = DiffRateA*t/Sc;		% Diffusion Rate A (D = L² / 2d*t)
Db = DiffRateB*t/Sc;		% Diffusion Rate B
Dr = Da/Db;					% Ratio of Da:Ds (1/Ls)^2;
Dn = Da/Dr;					% new D after scaling L
k = sqrt(dm*Da);			% stdev of D's step size distribution
L = sqrt(2*dm*Da);			% average diagonal (2D) step size
Lx = L/sqrt(2);             % average linear (1D) step size
Ls = 1/sqrt(Dr);			% scales Lx values for Dn
MSD = 2*dm*Da;				% mean squared displacement


XYL = zeros(2,Ndots);		% XY particle locations
XYS = zeros(2,Ndots);		% XY step sizes
XYLp = zeros(2,NSteps);		% preallocate matrix for trace dot

POLYSz = POLYSz./Scale;		% scale enclosures 
XWIDE = POLYSz(1)/2;		% half X enclosure size (will double below)
YHIGH = POLYSz(2)/2;		% half X enclosure size (will double below)

% [POLYSz(1)/2 POLYSz(2)/2 POLYSz(1) POLYSz(2)]
BOARDER = [-POLYSz(1)/2 -POLYSz(2)/2 POLYSz(1) POLYSz(2)];



XYL(1,:) = XYL(1,:)+5;
XYL(2,:) = XYL(2,:)+20;

% Box1 location [X Y W H]
BOXxy=10; BOXwh=14;
BOXLOC1 = [-POLYSz(1)/BOXxy POLYSz(2)/BOXxy*2 POLYSz(1)/BOXwh*2 POLYSz(2)/BOXwh*2]; 	
BOXLOC2 = [-POLYSz(1)/BOXxy -POLYSz(2)/BOXxy*4 POLYSz(1)/BOXwh*2 POLYSz(2)/BOXwh*2];

Box1LB = BOXLOC1(1:2);
Box1RT = Box1LB + BOXLOC1(3:4);
Box2LB = BOXLOC2(1:2);
Box2RT = Box2LB + BOXLOC2(3:4);


%================================================%
%               FIGURE SETUP
%------------------------------------------------%
Flh = figure(1);
set(Flh,'Units','pixels');  scnsize = get(0,'ScreenSize');
pos1 = [scnsize(3)/2  scnsize(4)/3  scnsize(3)/2  scnsize(4)/1.8];
set(Flh,'OuterPosition',pos1)
set(gcf,'Color',[.9,.9,.9])
%--------
xlim = [-XWIDE*1.01 XWIDE*1.01]; ylim = [-YHIGH*1.01 YHIGH*1.01];
%---
subplot('Position',[.05 .05 .65 .9])
Ph1 = scatter(XYL(1,:),XYL(2,:),5,[0 0 1]);
axis([xlim, ylim]);
% axis off
set(gca,'XTick',[],'YTick',[],'Color',[.8 .8 .8])
hold on

pause(1)
% SPh2 = subplot(1,1,1);
% Ph2 = plot(XYLp(1,:),XYLp(2,:));
% axis([xlim, ylim]);
% axis off
% hold(SPh2,'on');
rectangle('Position',BOXLOC1); hold on;
rectangle('Position',BOXLOC2); hold on;
rectangle('Position',BOARDER); hold on;

pause(1)

Ah3 = subplot('Position',[.75 .05 .22 .9]);
Ph3 = bar([1 1]);
set(Ah3,'Ylim',[0 50])

pause(1)
%-------------------------------------------------%



%===========================================%
for Nt = 1:NSteps 

	XYS = (k * randn(2,Ndots));	% generates step sizes
	
	% Sticky Functions
	%{.
	% [XYS] = STICKYFUN(XYL,XYS,Box1LB,Box2LB,Ls);
	Box1 = inboxfun(Box1LB,Box1RT,XYL);
	Box2 = inboxfun(Box2LB,Box2RT,XYL);
	XYS(:,Box1) = XYS(:,Box1)*(Ls);
	XYS(:,Box2) = XYS(:,Box2)*(Ls);
	%}
	
	XYL = XYL+XYS;				% adds step to location
	
	% Keep everything inside enclosure  %
	[XYL] = ENCLOSE(Nt,XYL,XWIDE,YHIGH,Ndots);

	
	doTrace=0;
	if doTrace
	XYLp(:,Nt) = XYL(:,1);		% save step of first dot (for trace)
	DUALPLOT(Nt,XYL,Ph1,XYLp,SPh2)
	XL(Nt,:) = XYL(1,:);
	YL(Nt,:) = XYL(2,:);
	end
	
	% figure(SPh2)
	set(Ph1,'XData',XYL(1,:),'YData',XYL(2,:));
	drawnow
	
	if mod(Nt,50)==0
	Box1N = numel(find(Box1>0));
	Box2N = numel(find(Box2>0));
	
	set(Ph3,'YData',[Box1N Box2N]);
	drawnow
	
	% if Nt>300;keyboard;end;
	
	end
	
%-------------------------------%
% pause(Delay)
if mod(Nt,100)==0;disp(Nt); end;
end % for Nt = 1:Nsteps 
%===========================================%


%===============================%
% varargout = {XYLp}; % export traced particle location
varargout = {XYL};
end
%===============================%



%-------------------------------------------%
% STICKYFUN: keep particles inside polygon
%-------------------------------------------%
function [XYS] = STICKYFUN(XYL,XYS,Box1LB,Box2LB,Ls)

INPSD1 = inboxfun(Box1LB,Box1RT,XYL);
INPSD2 = inboxfun(Box2LB,Box2RT,XYL);

XYS(:,INPSD1) = XYL(:,INPSD1)*(Ls);
XYS(:,INPSD2) = XYL(:,INPSD2)*(Ls);

end



%-------------------------------------------%
% ENCLOSE: keep particles inside polygon
%-------------------------------------------%
function [XYL] = ENCLOSE(Nt,XYL,XWIDE,YHIGH,Ndots)

	for j = 1:Ndots 
		if XYL(1,j)>(XWIDE) || XYL(1,j)<(-XWIDE)
			XYL(1,j) = sign(XYL(1,j))*(XWIDE);
		elseif XYL(2,j)>(YHIGH) || XYL(2,j)<(-YHIGH)
			XYL(2,j) = sign(XYL(2,j))*(YHIGH);
		end
	end

end



%===================================%
%		inboxfun
%===================================%
% Tests whether particles are in a box polygon
% and returns a logical vector
%===================================%
function [inbox] = inboxfun(LB,RT,xyl)

if LB(1)>RT(1)
	LBt=LB;
	RTt=RT;
	LB(1)=RTt(1);
	RT(1)=LBt(1);
end
if LB(2)>RT(2)
	LBt=LB;
	RTt=RT;
	LB(2)=RTt(2);
	RT(2)=LBt(2);
end

xylLB1 = xyl(1,:) > LB(1);
xylRT1 = xyl(1,:) < RT(1);
xylLB2 = xyl(2,:) > LB(2);
xylRT2 = xyl(2,:) < RT(2);

inbox = xylLB1 & xylRT1 & xylLB2 & xylRT2;

end
%===================================%



%-------------------------------------------%
% DUALPLOT: plots diffusion
%-------------------------------------------%
function [] = DUALPLOT(Nt,XYL,Ph1,XYLp,SPh2)
%-------------------------------------------%

LLg = 3;
if Nt > LLg
xp = (XYLp(1,(Nt-LLg):Nt));
yp = (XYLp(2,(Nt-LLg):Nt));
else
xp = (XYLp(1,Nt));
yp = (XYLp(2,Nt));
end

set(Ph1,'XData',XYL(1,:),'YData',XYL(2,:));
drawnow

plot(xp,yp,'Parent',SPh2);
drawnow

 
end



