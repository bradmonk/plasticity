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
Nsteps = 3000;	BRp = 50;
lotol = .2;		hitol = 1.9;
PSDsz = 8;		PSAsz = 4;
%---
vars = [Nsteps lotol hitol PSDsz PSAsz];
%-------------------------
doAMPARs = 1; AMPARN = 20; 
amparate = 1; G1RT = 1.1; G1ST = -3;
%---
glu = [doAMPARs AMPARN amparate G1RT G1ST];
%-------------------------

dT = .0014;

Lon = 2;
Bon = [0 10];
Ron = 15.0;

Loff = 2;
Boff = 1;
Roff = [0 10];


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

if dLon;	PLon  = linspace(Lon(1),Lon(2),BRp);
else		PLon = Lon.*ones(1,BRp);					end;
if dBon;	PBon  = linspace(Bon(1),Bon(2),BRp);
else		PBon = Bon.*ones(1,BRp);					end;
if dRon;	PRon  = linspace(Ron(1),Ron(2),BRp);
else		PRon = Ron.*ones(1,BRp);					end;
if dLoff;	PLoff = linspace(Loff(1),Loff(2),BRp);
else		PLoff = Loff.*ones(1,BRp);					end;
if dBoff;	PBoff = linspace(Boff(1),Boff(2),BRp);
else		PBoff = Boff.*ones(1,BRp);					end;
if dRoff;	PRoff = linspace(Roff(1),Roff(2),BRp);
else		PRoff = Roff.*ones(1,BRp);					end;



%==========================================================%
%==========================================================%
wbh = waitbar(0,'Initializing... ');
nn = 0;

%-------------------------
for aa = 1:BRp; 
	for bb = 1:BRp;
%-------------------------
	nn = nn+1;
	
	Lon = PLon(aa);
	Bon = PBon(aa);
	Ron = PRon(aa);
	Loff = PLoff(aa);
	Boff = PBoff(aa);
	Roff = PRoff(bb);


	Spara = [Lon Bon Ron Loff Boff Roff];
	[stepN OKGO] = RobustSAP4(Lon,Bon,Ron,Loff,Boff,Roff,dT,vars,glu);

	AOK(nn) = OKGO;
	LRPAR(nn,:) = Spara;
	STP(nn) = stepN;

	%-------------------------
	end;
	%-------------------------



%-------------------------
if mod(aa,10)==0;wbh = waitbar((aa/BRp),wbh,sprintf('%0.2f',(aa/BRp)*100));end;
end;
%==========================================================%
close(wbh)
%==========================================================%

keyboard


%====================================================%
%					DATA PREP
%====================================================%
%-----------
% LRPAR		: all the parameter combos that were run
%			: column headers: [Lon Bon Ron Loff Boff Roff]
%-----------
% didLBR	: 2 index for columns of the 2 test parameters 
%-----------

PMX = [LRPAR AOK' STP'];

didLBR = find(doLBR);
METASTEP = [STP' LRPAR(:,didLBR)];

STB8 = PMX(((AOK')==1),1:8);	% Stable P-combos
STB2 = STB8(:, didLBR);
DEC8 = PMX(((AOK')==0),1:8);	% Decay  P-combos
DEC2 = DEC8(:, didLBR);
EXP8 = PMX(((AOK')==2),1:8);	% Expand P-combos
EXP2 = EXP8(:, didLBR);


METAS = METASTEP;
METAS(((AOK')==2),1) = Nsteps+(Nsteps-METAS(((AOK')==2),1));



%====================================================%
%					FIGURES
%====================================================%
scsz = get(0,'ScreenSize');
xlabf1 = vnames(didLBR(2));	ylabf1 = vnames(didLBR(1));
xlabf2 = vnames(didLBR(2));	ylabf2 = vnames(didLBR(1));


%====================================================%
%			FIG1 - 2D Cluster Stability
%----------------------------------------------------%
fh1 = figure(1); hold on
set(fh1,'OuterPosition',(scsz./[2e-3 2e-2 2 1.5]))
%----------------------------------------------------%
scatter(STB2(:,2),STB2(:,1), sqrt(STB8(:,8)), 'ob','fill')
ylabel(ylabf1); xlabel(xlabf1);
hold on
scatter(DEC2(:,2),DEC2(:,1), sqrt(DEC8(:,8)), 'or','fill')
ylabel(ylabf1); xlabel(xlabf1);
hold on
scatter(EXP2(:,2),EXP2(:,1), sqrt(Nsteps*2), 'og','fill')
ylabel(ylabf1); xlabel(xlabf1);
title('Stable(blue) - Decay(red) - Expand(green)')
%----------------------------------------------------%



%====================================================%
%			FIG2 - 3D Cluster Stability
%----------------------------------------------------%
fh2 = figure(2); hold on
set(fh2,'OuterPosition',(scsz./[2e-3 2e-3 2 1.5]))
%----------------------------------------------------%
zz = METASTEP(:,1);
xx = METASTEP(:,3);
yy = METASTEP(:,2);

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
grid on
colormap jet
view(-15, 65);
hXLabel = xlabel(xlabf2);
hYLabel = ylabel(ylabf2);
hZLabel = zlabel('stepN');
set(gca,'FontName','Helvetica');
set([hXLabel, hYLabel, hZLabel],'FontName','AvantGarde');
set([hXLabel, hYLabel, hZLabel],'FontSize',12);
%----------------------------------------------------%


%====================================================%
%			FIG3 - 3D Cluster Life
%----------------------------------------------------%
fh3 = figure(3); hold on
set(fh3,'OuterPosition',(scsz./[2e-3 2e-3 2 1.5]))
%----------------------------------------------------%
zz = METAS(:,1);
xx = METAS(:,3);
yy = METAS(:,2);

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
grid on
colormap jet
view(140, 40);
hXLabel = xlabel(xlabf2);
hYLabel = ylabel(ylabf2);
hZLabel = zlabel('stepN');
set(gca,'FontName','Helvetica');
set([hXLabel, hYLabel, hZLabel],'FontName','AvantGarde');
set([hXLabel, hYLabel, hZLabel],'FontSize',12);
%----------------------------------------------------%

%====================================================%
%			FIG11 - 3D Cluster Stability
%----------------------------------------------------%
%{
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


METASTABLE = LRPAR((AOK==1),1:6);
DECAY = LRPAR((AOK==0),1:6);
EXPAND = LRPAR((AOK==2),1:6);
METASTB = METASTABLE(:, didLBR);	% Stable P-combos
DEC = DECAY(:, didLBR);				% Decay  P-combos
EXP = EXPAND(:, didLBR);			% Expand P-combos


%====================================================%
%			FIG1 - 2D Cluster Stability
%----------------------------------------------------%
fh1 = figure(1); hold on
set(fh1,'OuterPosition',(scsz./[2e-3 2e-2 2 1.5]))
%----------------------------------------------------%
scatter(STB2(:,2),STB2(:,1), sqrt(STB8(:,8)), 'ob','fill')
ylabel(ylabf1); xlabel(xlabf1);
hold on
scatter(DEC2(:,2),DEC2(:,1), sqrt(DEC8(:,8)), 'or','fill')
ylabel(ylabf1); xlabel(xlabf1);
hold on
scatter(EXP2(:,2),EXP2(:,1), sqrt((Nsteps-EXP8(:,8)).*2+Nsteps), 'og','fill')
ylabel(ylabf1); xlabel(xlabf1);
title('Stable(blue) - Decay(red) - Expand(green)')
%----------------------------------------------------%

%====================================================%
%			FIG1 - 2D Cluster Stability
%----------------------------------------------------%
fh1 = figure(1); hold on
set(fh1,'OuterPosition',(scsz./[2e-3 2e-2 2 1.5]))
%----------------------------------------------------%
scatter(METASTB(:,2),METASTB(:,1), (30), 'ob','fill')
ylabel(ylabf1); xlabel(xlabf1);
hold on
scatter(DEC(:,2),DEC(:,1), (10), 'or','fill')
ylabel(ylabf1); xlabel(xlabf1);
hold on
scatter(EXP(:,2),EXP(:,1), (50), 'og','fill')
ylabel(ylabf1); xlabel(xlabf1);
title('Stable(blue) - Decay(red) - Expand(green)')
%----------------------------------------------------%


%====================================================%
%			FIG12 - 3D Cluster Stability
%----------------------------------------------------%
fig12 = figure(12); set(12,'Units','pixels'); scsz = get(0,'ScreenSize');
set(fig12,'OuterPosition',[scsz(3)/3  scsz(4)/5  scsz(3)/2  scsz(4)/1.5]);
xlab = vnames(didLBR(2));
ylab = vnames(didLBR(1));
%----------------------------------------------------%
figure(12); hold on
zz = METASTEP(:,1);
xx = METASTEP(:,3);
yy = METASTEP(:,2);


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
view(-15, 65);
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
%}




