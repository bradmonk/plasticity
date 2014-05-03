clc, close all, clear all

% QRpoints: number of points to test in range specified below
QRpoints = 5;
NPtested = 3;
assignin('base','QRpoints',QRpoints);

% % RATE CONSTANTS
% Lon = 2.0;
% Loff = 2.0;
% %-- INDY (low=fast)
% %Qon = [10 15];
% Qon = [7 7];
% Qoff = [.1 5];
% %-- DEP (high=fast)
% Ron = [.1 10];
% Roff = [.1 3];


% RATE CONSTANTS
Lon = 2.0;
Loff = 2.0;
%-- INDY (low=fast)
%Qon = [10 15];
Qon = [7 15];
Qoff = [1 5];
%-- DEP (high=fast)
Ron = [.1 10];
Roff = [.1 2.5];


% profile on
%==========================================================%
% Sfull{25,50} = [];
% Lon = 1.8;  % On Energy (lower = more on events)
% Qon = 32;   % On Neighbor-Independant Rate (new growth) (lower = more on)
% Ron = 3;    % On Neighbor-Dependant Rate (cluster fill-in) (higher = more on)
%  
% Loff = 1.8; % Off Energy (higher = more off events)
% Qoff = 1.2; % Off Neighbor-Inependant Rate (uniform off)  (lower = more off)
% Roff = 3;   % Off Neighbor-Dependant Rate (edge off) (higher = more off)
%==========================================================%

PQon  = linspace(Qon(1),Qon(2),QRpoints);
PQoff = linspace(Qoff(1),Qoff(2),QRpoints);
PRon  = linspace(Ron(1),Ron(2),QRpoints);
PRoff = linspace(Roff(1),Roff(2),QRpoints);



nn = 0;
% for aa = 1:QRpoints; 
for bb = 1:QRpoints; for cc = 1:QRpoints; for dd = 1:QRpoints;
%disp(aa);
disp(bb);disp(cc);disp(dd);
nn = nn+1;

	% Qon = PQon(aa);
	Qon = PQon(1);
	Ron = PQoff(bb);

	Qoff = PRon(cc);
	Roff = PRoff(dd);

	Spara = [Lon Loff Qon Qoff Ron Roff];
	[Scell HKcell OKGO] = SAPbrute(Lon,Qon,Ron,Loff,Qoff,Roff);

	AOK(nn) = OKGO;
	LRPAR(nn,:) = Spara;
	Sfull{nn} = {Spara Scell OKGO};
	HKfull{nn} = {Spara HKcell OKGO};
	
%{	
figure(1);
subplot(2,2,1),imagesc(Scell{1});
subplot(2,2,2),imagesc(Scell{2});
subplot(2,2,3),imagesc(Scell{3});
subplot(2,2,4),imagesc(Scell{4});
colormap('bone')
title('2D Particle Map');
title('PSD1');
drawnow;
%}
%==========================================================%
end;end;end;
%end;
%==========================================================%
% profile viewer
%==========================================================%
save('RobustData.mat', 'Sfull','HKfull')
% clear
% load RobustData.mat
%==========================================================%

OKPARAMS = LRPAR((AOK>0),3:6);
NOKPARAMS = LRPAR((AOK<1),3:6);

numel(OKPARAMS(:,1))

min(OKPARAMS)
max(OKPARAMS)

%=============================%
%-----------------------------%
figure(10); hold on
X = [OKPARAMS(:,2) OKPARAMS(:,3)];
V = OKPARAMS(:,4);
plot3(X(:,1),X(:,2),(OKPARAMS(:,1).*0), '.b')
grid
stem3(X(:,1),X(:,2),V,'sb','fill')
view(-36, 60);

X = [NOKPARAMS(:,2) NOKPARAMS(:,3)];
V = NOKPARAMS(:,4);
plot3(X(:,1),X(:,2),(NOKPARAMS(:,1).*0), '.r')
stem3(X(:,1),X(:,2),V,':sr')
view(-36, 60);
%=============================%
figure(11)
plot3(OKPARAMS(:,2),OKPARAMS(:,3),OKPARAMS(:,4), '.b')
grid; hold on
plot3(NOKPARAMS(:,2),NOKPARAMS(:,3),NOKPARAMS(:,4), '.r')
%=============================%


%=============================%

%{
[PQonX,PQonY] = meshgrid(PQon);
[PQoffX,PQoffY] = meshgrid(PQoff);
[PRonX,PRonY] = meshgrid(PRon);
[PRoffX,PRoffY] = meshgrid(PRoff);


TOEPMx = toeplitz(PQoff,PRon);
mesh(PQonX,PQonY,TOEPMx)

% SLICE 3D VOLUME
[x,y,z] = meshgrid(PQon,PQoff,PRon);
FBFBT = ones(4,4,4);
FBFBT(:,:,1) = (PRoffX.*FBFBT(:,:,1));
FBFBT(:,:,2) = (PRoffX.*FBFBT(:,:,2));
FBFBT(:,:,3) = (PRoffX.*FBFBT(:,:,3));
FBFBT(:,:,4) = (PRoffX.*FBFBT(:,:,4));
v = FBFBT;
xslice = [PQon(2),PQon(3)]; yslice = PQoff(2); zslice = [PRon(2),PRon(3)];
slice(x,y,z,v,xslice,yslice,zslice)
colormap hsv



%=============================%
[x,y,z] = meshgrid(-2:.2:2,-2:.25:2,-2:.16:2);
v = x.*exp(-x.^2-y.^2-z.^2);
xslice = [-1.2,.8,2]; yslice = 2; zslice = [-2,0];
slice(x,y,z,v,xslice,yslice,zslice)
colormap hsv
for i = -2:.5:2
 hsp = surf(linspace(-2,2,20),linspace(-2,2,20),zeros(20)+i);
 rotate(hsp,[1,-1,1],30)
 xd = get(hsp,'XData');
 yd = get(hsp,'YData');
 zd = get(hsp,'ZData');
 delete(hsp)
 slice(x,y,z,v,[-2,2],2,-2) % Draw some volume boundaries
 hold on
 slice(x,y,z,v,xd,yd,zd)
 hold off
 axis tight 
 view(-5,10) 
 drawnow
 pause(.5)
end
%-----------------------------%
[xsp,ysp,zsp] = sphere;
slice(x,y,z,v,[-2,2],2,-2)  % Draw some volume boundaries
for i = -3:.2:3
 hsp = surface(xsp+i,ysp,zsp);
 rotate(hsp,[1 0 0],90)
 xd = get(hsp,'XData');
 yd = get(hsp,'YData');
 zd = get(hsp,'ZData');
 delete(hsp)
 hold on
 hslicer = slice(x,y,z,v,xd,yd,zd);
 axis tight 
 xlim([-3,3])
 view(-10,35)
 drawnow
 pause(.5)
 delete(hslicer)
 hold off
end
%-----------------------------%
PQonScat = [PQon' PQoff';PQon' PRon';PQon' PRoff';...
	PQoff' PRon'; PQoff' PRoff'; PRon' PRoff'];
X = PQonScat;
V = PQonScat(:,1);
hold on
plot3(X(:,1),X(:,2),zeros(numel(PQonScat(:,1)),1), '*r')
grid
stem3(X(:,1),X(:,2),V,'^','fill')
hold off
view(322.5, 30);




%=============================%
% Cell{ID}{data}{step}
%-----------------------------%
% ID = 1:QRpoints^QRpoints
%---
% data = 1 = param vals
% data = 2 = SAP Mx
% data = 3 = Robust yes/no indicator
%---
% step = 1 = step 10
% step = 2 = step 100
% step = 3 = step 1000
% step = 4 = step 5000
%---
% cellplot(Sfull{1})
%=============================%

Sfull{1}{3}




sum(Sfull{1}{2}{1})
sum(HKfull{1}{2}{1})

imagesc(Sfull{1}{2}{4})
imagesc(HKfull{1}{2}{4})

Sf = Sfull{1}{2}{4};
Sf1 = numel(find(Sf(1,1:end)));
Sf2 = numel(find(Sf(1:end,1)));

Sfull{1}{2}{4}
HKfull{1}{2}{4}

%}




% cellplot(Sfull{1})
%==========================================================%
doMakeMovie = 0;
if doMakeMovie
%==========================================================%
fig1 = figure(1);
set(gcf,'OuterPosition',[800,300,500,500])

for kk = 1:QRpoints^NPtested
figure(1);
subplot('Position',[.05 .55 .40 .40]),imagesc(Sfull{kk}{2}{1});
set(gca,'XTick',[],'YTick',[])
subplot('Position',[.55 .55 .40 .40]),imagesc(Sfull{kk}{2}{2});
set(gca,'XTick',[],'YTick',[])
subplot('Position',[.05 .05 .40 .40]),imagesc(Sfull{kk}{2}{3});
set(gca,'XTick',[],'YTick',[])
subplot('Position',[.55 .05 .40 .40]),imagesc(Sfull{kk}{2}{4});
set(gca,'XTick',[],'YTick',[])
colormap('bone')
text(1,-2,['Lon       Loff       Qon       Qoff       Ron       Roff'],...
	'FontSize',16,'HorizontalAlignment','center')
text(1,-.5,['' num2str(Sfull{kk}{1}) ''],...
	'FontSize',16,'HorizontalAlignment','center')
Vid(kk) = getframe(gcf,[0 0 500 500]);
% writeVideo(ClusVid,Vid);
end
%==========================================================%
% PLAY MOVIE
Vidframes = numel(Vid);

if Vidframes >= 150
fig2 = figure(2);
set(gcf,'OuterPosition',[800,300,500,500])
movie(fig2,Vid(1:150),1,5)
movie(fig2,Vid(151:end),1,5)
else
fig2 = figure(2);
set(gcf,'OuterPosition',[800,300,500,500])
movie(fig2,Vid(1:end),1,5)
end

%==========================================================%
end % if doMakeMovie
%==========================================================%
% %==========================================================%
% doMakeMovie = 0;
% if doMakeMovie
% %==========================================================%
% fig1 = figure(1);
% set(gcf,'OuterPosition',[800,300,500,500])
% 
% for kk = 1:NPtested^QRpoints
% figure(1);
% subplot('Position',[.05 .55 .40 .40]),imagesc(HKfull{kk}{2}{1});
% set(gca,'XTick',[],'YTick',[])
% subplot('Position',[.55 .55 .40 .40]),imagesc(HKfull{kk}{2}{2});
% set(gca,'XTick',[],'YTick',[])
% subplot('Position',[.05 .05 .40 .40]),imagesc(HKfull{kk}{2}{3});
% set(gca,'XTick',[],'YTick',[])
% subplot('Position',[.55 .05 .40 .40]),imagesc(HKfull{kk}{2}{4});
% set(gca,'XTick',[],'YTick',[])
% colormap('bone')
% text(1,-2,['Lon       Loff       Qon       Qoff       Ron       Roff'],...
% 	'FontSize',16,'HorizontalAlignment','center')
% text(1,-.5,['' num2str(Sfull{kk}{1}) ''],...
% 	'FontSize',16,'HorizontalAlignment','center')
% Vid(kk) = getframe(gcf,[0 0 500 500]);
% % writeVideo(ClusVid,Vid);
% end
% %==========================================================%
% % PLAY MOVIE
% fig2 = figure(2);
% set(gcf,'OuterPosition',[800,300,500,500])
% movie(fig2,Vid(1:150),1,15)
% movie(fig2,Vid(151:end),1,15)
% %==========================================================%
% end % if doMakeMovie
% %==========================================================%



