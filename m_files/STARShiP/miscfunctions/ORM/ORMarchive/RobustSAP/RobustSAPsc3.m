clc, close all, clear all



%{
% Lon = [0.5 3.0];
% Bon = [0.1 80.0];
% Ron = 9.5;
% 
% Loff = 2.5;
% Boff = 0;
% Roff = 2;


Lon = [ 0.5  3.5  ];
Bon = [ 0.1  80.0 ];
Ron = 15;

Loff = 2;
Boff = 1;
Roff = 4;
%}
%==========================================================%
BRpoints = 40;	NPtested = 3;
lotol = .2;		hitol = 1.5;
Nsteps = 2000;	doAMPARs = 0;
PSDsz = 8;		PSAsz = 4;
%---
vars = [Nsteps lotol hitol doAMPARs PSDsz PSAsz];
%-------------------------

dT = .0014;

Lon = 2;
Bon = [0 10];
Ron = 15.0;

Loff = 2;
Boff = 1;
Roff = [1.0 15.0];


%--		[Lon	Bon		Ron		Loff	Boff	Roff]	--%
doLBR = [0		1		0		0		0		1	];

%--------------
varlabs = {'dLon'; 'dBon'; 'dRon'; 'dLoff'; 'dBoff'; 'dRoff'};
varsl = {varlabs};

vnames = genvarname(varsl{1}, who);

for m=1:numel(doLBR)
	eval([vnames{m} ' = doLBR(m);']); 
end
%--------------

if dLon;	PLon  = linspace(Lon(1),Lon(2),BRpoints);
else		PLon = Lon.*ones(1,BRpoints);					end;
if dBon;	PBon  = linspace(Bon(1),Bon(2),BRpoints);
else		PBon = Bon.*ones(1,BRpoints);					end;
if dRon;	PRon  = linspace(Ron(1),Ron(2),BRpoints);
else		PRon = Ron.*ones(1,BRpoints);					end;
if dLoff;	PLoff = linspace(Loff(1),Loff(2),BRpoints);
else		PLoff = Loff.*ones(1,BRpoints);					end;
if dBoff;	PBoff = linspace(Boff(1),Boff(2),BRpoints);
else		PBoff = Boff.*ones(1,BRpoints);					end;
if dRoff;	PRoff = linspace(Roff(1),Roff(2),BRpoints);
else		PRoff = Roff.*ones(1,BRpoints);					end;





%-------------------------
wbh = waitbar(0,'Initializing... ');
nn = 0;
for aa = 1:BRpoints; 
for bb = 1:BRpoints; 
nn = nn+1;
%-------------------------

	Lon = PLon(aa);
	Bon = PBon(aa);
	Ron = PRon(aa);
	Loff = PLoff(aa);
	Boff = PBoff(aa);
	Roff = PRoff(bb);


	Spara = [Lon Bon Ron Loff Boff Roff];
	[stepN OKGO] = RobustSAP2(Lon,Bon,Ron,Loff,Boff,Roff,dT,vars);

	AOK(nn) = OKGO;
	LRPAR(nn,:) = Spara;
	STP(nn) = stepN;

%-------------------------
end;
%-------------------------

	if mod(aa,10)==0
    perc = aa/BRpoints;
    wbh = waitbar(perc,wbh,sprintf('%0.2f',perc*100));
	end

%-------------------------
end;
%==========================================================%
close(wbh)
%==========================================================%



%-- [Lon Loff Bon Boff Ron Roff] --%
METASTABLE = LRPAR((AOK==1),1:6);
DECAY = LRPAR((AOK==0),1:6);
EXPAND = LRPAR((AOK==2),1:6);

%--	doLBR = [Lon Bon Ron Loff Boff Roff] --%
didLBR = find(doLBR);

METASTB = METASTABLE(:, didLBR);
DEC = DECAY(:, didLBR);
EXP = EXPAND(:, didLBR);

xlab = vnames(didLBR(1));
ylab = vnames(didLBR(2));

%====================================================%
%				FIGURE 10 SETUP
%----------------------------------------------------%
fig10 = figure(10); set(10,'Units','pixels'); scsz = get(0,'ScreenSize');
set(fig10,'OuterPosition',[scsz(3)/3  scsz(4)/5  scsz(3)/2.5  scsz(4)/2]);
%----------------------------------------------------%
figure(10); hold on
scatter(METASTB(:,2),METASTB(:,1), 'og','fill')
ylabel(ylab); xlabel(xlab);
hold on
scatter(DEC(:,2),DEC(:,1), 'ob','fill')
ylabel(ylab); xlabel(xlab);
hold on
scatter(EXP(:,2),EXP(:,1), 'or','fill')
ylabel(ylab); xlabel(xlab);
%----------------------------------------------------%

%-- LRPAR = [Lon Bon Ron Loff Boff Roff] --%
LBRO = LRPAR(:,didLBR);
METASTEP = [STP' LBRO(:,1) LBRO(:,2)]; 

%====================================================%
%				FIGURE 11 SETUP
%----------------------------------------------------%
fig11 = figure(11); set(11,'Units','pixels'); scsz = get(0,'ScreenSize');
set(fig11,'OuterPosition',[scsz(3)/3  scsz(4)/5  scsz(3)/2  scsz(4)/1.5]);
%----------------------------------------------------%
figure(11); hold on
X = [METASTEP(:,3) METASTEP(:,2)];
V = METASTEP(:,1);
plot3(X(:,1),X(:,2),(METASTEP(:,1)), '.b')
% stem3(X(:,1),X(:,2),V,':ob','fill')
grid on
view(-12, 79);
hXLabel = xlabel(xlab);
hYLabel = ylabel(ylab);
hZLabel = zlabel('stepN');
set(gca,'FontName','Helvetica');
set([hXLabel, hYLabel, hZLabel],'FontName','AvantGarde');
set([hXLabel, hYLabel, hZLabel],'FontSize',12);
%----------------------------------------------------%




%====================================================%
%				FIGURE 12 SETUP
%----------------------------------------------------%
fig12 = figure(12); set(12,'Units','pixels'); scsz = get(0,'ScreenSize');
set(fig12,'OuterPosition',[scsz(3)/3  scsz(4)/5  scsz(3)/2  scsz(4)/1.5]);
%----------------------------------------------------%
figure(12); hold on
zz = METASTEP(:,1);
xx = METASTEP(:,2);
yy = METASTEP(:,3);


F = TriScatteredInterp(xx,yy,zz);

xxmin = min(xx);
yymin = min(yy);
xxmax = max(xx);
yymax = max(yy);

til = yymin:((yymax-yymin)/80):yymax;
tih = xxmin:((xxmax-xxmin)/80):xxmax;
[qx,qy] = meshgrid(tih,til);
qz = F(qx,qy);
mesh(qx,qy,qz);
% surf(qx,qy,qz);
grid on
colormap jet
view(-15, 75);
hXLabel = xlabel(xlab);
hYLabel = ylabel(ylab);
hZLabel = zlabel('stepN');
set(gca,'FontName','Helvetica');
set([hXLabel, hYLabel, hZLabel],'FontName','AvantGarde');
set([hXLabel, hYLabel, hZLabel],'FontSize',12);
% hold on;
% plot3(xx,yy,zz,'ro');
% view(-15, 75);
%----------------------------------------------------%





%{
%====================================================%
%				FIGURE 10 SETUP
%----------------------------------------------------%
fig10 = figure(10); set(10,'Units','pixels'); scsz = get(0,'ScreenSize');
set(fig10,'OuterPosition',[scsz(3)/3  scsz(4)/5  scsz(3)/2.5  scsz(4)/2]);
%----------------------------------------------------%
figure(10); hold on
scatter(METASTABLE(:,2),METASTABLE(:,1), 'og','fill')
ylabel('Lon'); xlabel('Bon');
hold on
scatter(DECAY(:,2),DECAY(:,1), 'ob','fill')
ylabel('Lon'); xlabel('Bon');
hold on
scatter(EXPAND(:,2),EXPAND(:,1), 'or','fill')
ylabel('Lon'); xlabel('Bon');
%----------------------------------------------------%



%-- LRPAR = [Lon Bon Ron Loff Boff Roff] --%
%-- METASTEP = [stepN Bon Lon] --%

METASTEP = [STP' LRPAR(:,1) LRPAR(:,2)]; 


%====================================================%
%				FIGURE 11 SETUP
%----------------------------------------------------%
fig11 = figure(11); set(11,'Units','pixels'); scsz = get(0,'ScreenSize');
set(fig11,'OuterPosition',[scsz(3)/3  scsz(4)/5  scsz(3)/2  scsz(4)/1.5]);
%----------------------------------------------------%
figure(11); hold on
X = [METASTEP(:,3) METASTEP(:,2)];
V = METASTEP(:,1);
plot3(X(:,1),X(:,2),(METASTEP(:,1)), '.b')
% stem3(X(:,1),X(:,2),V,':ob','fill')
grid on
view(-12, 79);
hXLabel = xlabel('Lon');
hYLabel = ylabel('Bon');
hZLabel = zlabel('stepN');
set(gca,'FontName','Helvetica');
set([hXLabel, hYLabel, hZLabel],'FontName','AvantGarde');
set([hXLabel, hYLabel, hZLabel],'FontSize',12);
%----------------------------------------------------%




%====================================================%
%				FIGURE 12 SETUP
%----------------------------------------------------%
fig12 = figure(12); set(12,'Units','pixels'); scsz = get(0,'ScreenSize');
set(fig12,'OuterPosition',[scsz(3)/3  scsz(4)/5  scsz(3)/2  scsz(4)/1.5]);
%----------------------------------------------------%
figure(12); hold on
xx = METASTEP(:,3);
yy = METASTEP(:,2);
zz = METASTEP(:,1);

F = TriScatteredInterp(xx,yy,zz);

xxmin = min(xx);
yymin = min(yy);
xxmax = max(xx);
yymax = max(yy);

til = yymin:((yymax-yymin)/80):yymax;
tih = xxmin:((xxmax-xxmin)/80):xxmax;
[qx,qy] = meshgrid(tih,til);
qz = F(qx,qy);
mesh(qx,qy,qz);
% surf(qx,qy,qz);
colormap jet
view(-15, 75);
hXLabel = xlabel('Lon');
hYLabel = ylabel('Bon');
hZLabel = zlabel('stepN');
set(gca,'FontName','Helvetica');
set([hXLabel, hYLabel, hZLabel],'FontName','AvantGarde');
set([hXLabel, hYLabel, hZLabel],'FontSize',12);
% hold on;
% plot3(xx,yy,zz,'ro');
% view(-15, 75);
%----------------------------------------------------%

%}






%%
% cellplot(Sfull{1})
%==========================================================%
doMakeMovie = 0;
if doMakeMovie
%==========================================================%
fig1 = figure(1);
set(gcf,'OuterPosition',[800,300,500,500])

for kk = 1:BRpoints^NPtested
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

















