clc, close all, clear all

Test2 = 1;
Test3 = 0;



%==========================================================%
if Test2
%-------------------------
QRpoints = 50;	NPtested = 3;

% RATE CONSTANTS
%-- B: INDY (low=fast) R: DEP (high=fast)

Lon = [0.5 3.0];
Bon = [0.1 80.0];
Ron = 9.5;

Loff = 2.5;
Boff = 0;
Roff = 2;

PLon  = linspace(Lon(1),Lon(2),QRpoints);
PBon  = linspace(Bon(1),Bon(2),QRpoints);

% PLoff = linspace(Loff(1),Loff(2),QRpoints);
% PBoff = linspace(Boff(1),Boff(2),QRpoints);
% PRon  = linspace(Ron(1),Ron(2),QRpoints);
% PRoff = linspace(Roff(1),Roff(2),QRpoints);
%-------------------------
h = waitbar(0,'Initializing waitbar...');
nn = 0;
for aa = 1:QRpoints; 
for bb = 1:QRpoints; 
% disp(aa);disp(bb);
nn = nn+1;
%-------------------------

	Lon = PLon(aa);
	Bon = PBon(bb);


	Spara = [Lon Loff Bon Boff Ron Roff];
	[Scell HKcell OKGO] = SRobustSAP(Lon,Bon,Ron,Loff,Boff,Roff);

	AOK(nn) = OKGO;
	LRPAR(nn,:) = Spara;
	Sfull{nn} = {Spara Scell OKGO};
	HKfull{nn} = {Spara HKcell OKGO};

%-------------------------
end;
%-------------------------

	if mod(aa,5)==0
    perc = aa/QRpoints;
    waitbar(perc,h,sprintf('%0.2f',perc*100))
	end

%-------------------------
end;
%==========================================================%
save('RobustData.mat', 'Sfull','HKfull')
end % if Test2
%==========================================================%


%-- [Lon Loff Bon Boff Ron Roff] --%
METASTABLE = LRPAR((AOK==1),1:6);
DECAY = LRPAR((AOK==0),1:6);
EXPAND = LRPAR((AOK==2),1:6);



%====================================================%
%							FIGURE SETUP
%----------------------------------------------------%
fig10 = figure(10); set(10,'Units','pixels'); scsz = get(0,'ScreenSize');
set(fig10,'OuterPosition',[scsz(3)/3  scsz(4)/5  scsz(3)/2.5  scsz(4)/2]);
%----------------------------------------------------%
figure(10); hold on
scatter(METASTABLE(:,3),METASTABLE(:,1), '.g')
ylabel('Lon'); xlabel('Bon');
hold on
scatter(DECAY(:,3),DECAY(:,1), '.b')
ylabel('Lon'); xlabel('Bon');
hold on
scatter(EXPAND(:,3),EXPAND(:,1), '.r')
ylabel('Lon'); xlabel('Bon');
%----------------------------------------------------%














%==========================================================%
if Test3
%-------------------------
NPtested = 3;

% RATE CONSTANTS
%-- B: INDY (low=fast) R: DEP (high=fast)
Lon = [1.8 2.4];	Loff = [2.5 2.9];
Bon = [1 20];		Boff = [0 0];
Ron = [9.5 9.5];	Roff = [2 2];


PLon  = linspace(Lon(1),Lon(2),QRpoints);
PLoff = linspace(Loff(1),Loff(2),QRpoints);

PBon  = linspace(Bon(1),Bon(2),QRpoints);
PBoff = linspace(Boff(1),Boff(2),QRpoints);

PRon  = linspace(Ron(1),Ron(2),QRpoints);
PRoff = linspace(Roff(1),Roff(2),QRpoints);
%-------------------------
nn = 0;
for aa = 1:QRpoints; 
for bb = 1:QRpoints; 
for cc = 1:QRpoints;
disp(aa);disp(bb);disp(cc);
nn = nn+1;
%-------------------------

	Bon = PBon(aa);
	Boff = PBoff(bb);
	
	Ron = PRon(cc);
	Roff = PRoff(cc);
	
	Lon = PLon(cc);
	Loff = PLoff(cc);


	Spara = [Lon Loff Bon Boff Ron Roff];
	[Scell HKcell OKGO] = SRobustSAP(Lon,Bon,Ron,Loff,Boff,Roff);

	AOK(nn) = OKGO;
	LRPAR(nn,:) = Spara;
	Sfull{nn} = {Spara Scell OKGO};
	HKfull{nn} = {Spara HKcell OKGO};

%-------------------------
end;end;end;
%==========================================================%
save('RobustData.mat', 'Sfull','HKfull')
end % if Test3
%==========================================================%







%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Test3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%==========================================================%
assignin('base','QRpoints',QRpoints);
% clear; load RobustData.mat;
%==========================================================%

%-- [Bon Boff Ron Roff PPRrat] --%
OKPARAMS = LRPAR((AOK>0),1:6);
NOKPARAMS = LRPAR((AOK<1),1:6);

%======================================================================%
%							FIGURE SETUP
%----------------------------------------------------------------------%
%----------------------------------------------------------------------%
fig10 = figure(10);
set(10,'Units','pixels');scnsize = get(0,'ScreenSize');
pos1 = [scnsize(3)/3  scnsize(4)/5  scnsize(3)/2.5  scnsize(4)/2];
set(fig10,'OuterPosition',pos1)
set(gcf,'Color',[.9,.9,.9])
%----------------------------------------------------------------------%

figure(10); hold on




X = [OKPARAMS(:,1) OKPARAMS(:,3)];
V = OKPARAMS(:,2);
% plot3(X(:,1),X(:,2),(OKPARAMS(:,1).*0), '.b')
plot3(X(:,1),X(:,2),(OKPARAMS(:,2)), '.b')
grid
stem3(X(:,1),X(:,2),V,':ob','fill')

X = [NOKPARAMS(:,1) NOKPARAMS(:,3)];
V = NOKPARAMS(:,2);
% plot3(X(:,1),X(:,2),(NOKPARAMS(:,1).*0), '.r')
%plot3(X(:,1),X(:,2),(NOKPARAMS(:,5)), '.r')
stem3(X(:,1),X(:,2),V,':or','fill')
grid on
view(-10, 11);
hXLabel = xlabel('Lon');
hYLabel = ylabel('Bon');
hZLabel = zlabel(num2str(['Loff (@ ' num2str(Loff) ')']) );
set(gca,'FontName','Helvetica');
set([hXLabel, hYLabel, hZLabel],'FontName','AvantGarde');
set([hXLabel, hYLabel, hZLabel],'FontSize',12);
%=============================%
% figure(10); hold on
% plot3(OKPARAMS(:,1),OKPARAMS(:,2),OKPARAMS(:,5), '.b')
% grid; hold on
% plot3(NOKPARAMS(:,1),NOKPARAMS(:,2),NOKPARAMS(:,5), '.r')
% view(-36, 60);
%=============================%

%=============================%
x = OKPARAMS(:,1);
y = OKPARAMS(:,3);
z = OKPARAMS(:,2);
xyz = [x y z];

% plot3(x,y,z,'.-')
tri = delaunay(x,y);
% plot(x,y,'.')
[r,c] = size(tri);
disp(r)
%=============================%
tri3 = delaunay(x,y,z);
figure(fig10)
h = trisurf(tri3, x, y, z);
axis vis3d
l = light('Position',[-50 -15 29]);
lighting phong
shading interp
hXLabel = xlabel('Lon');
hYLabel = ylabel('Bon');
hZLabel = zlabel(num2str(['Loff (@ ' num2str(Loff) ')']) );
%=============================%
% lighting flat
% lighting gouraud
% lighting phong
% lighting none
% shading flat 
% shading faceted 
% shading interp 
% Create light
%======================================================================%
set(gcf, 'PaperPositionMode', 'auto');
saveas(gcf,'SPACEPLOT','png');
% saveas(gcf, ['outputfigs/SPACEPLOT.png']);
%======================================================================%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end %if Test3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%









%%
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
text(1,-2,['Lon       Loff       Bon       Boff       Ron       Roff'],...
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

















