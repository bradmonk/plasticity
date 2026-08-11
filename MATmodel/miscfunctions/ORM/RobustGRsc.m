clc, close all, clear all

%==========================================================%
Nsteps = 3600;	BRp = 40;
lotol = .5;		hitol = 1.5;
PSDsz = 15;		PSAsz = 5;
doGaussianMask = 0;
%---
vars = [Nsteps lotol hitol PSDsz PSAsz doGaussianMask];
%-------------------------
doAMPARs = 0; AMPARN = 20; amparate = 1;
glu = [doAMPARs AMPARN amparate];
%-------------------------

dT = .0014;

Lon = 2;
Bon = 14;
Ron = [10 30];

Loff = 2;
Boff = 1;
Roff = [2 10];


GTon = 0;
GToff = -0;


%--		[Lon	Bon		Ron		Loff	Boff	Roff	GTon	GToff]	--%
doLBRG = [0		0		1		0		0		1		0		0	 ];

%--------------
varlabs = {'dLon';'dBon';'dRon';'dLoff';'dBoff';'dRoff';'dGTon';'dGToff'};
varsl = {varlabs};

vnames = genvarname(varsl{1}, who);

for m=1:numel(doLBRG)
	eval([vnames{m} ' = doLBRG(m);']); 
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
if dGTon;	PGTon  = linspace(GTon(1),GTon(2),BRp);
else		PGTon = GTon.*ones(1,BRp);					end;
if dGToff;	PGToff  = linspace(GToff(2),GToff(1),BRp);
else		PGToff = GToff.*ones(1,BRp);				end;


%==========================================================%
% letLBR = struct('dLon',{{'aa'}},'dBon',{{'aa'}},'dRon',{{'aa'}},...
% 	'dLoff',{{'aa'}},'dBoff',{{'aa'}},'dRoff',{{'aa'}},...
% 	'dGTon',{{'aa'}},'dGToff',{{'aa'}});
% letLBR.dLon{:}

numLBRG = struct('dLon',{PLon},'dBon',{PBon},'dRon',{PRon},...
	'dLoff',{PLoff},'dBoff',{PBoff},'dRoff',{PRoff},...
	'dGTon',{PGTon},'dGToff',{PGToff});

%==========================================================%
wbh = waitbar(0,'Initializing... ');
nn = 0;

%-------------------------
for aa = 1:BRp; 
	for bb = 1:BRp;
%-------------------------
	nn = nn+1;
	
	Lon = numLBRG.dLon(aa);
	Bon = numLBRG.dBon(aa);
	Ron = numLBRG.dRon(aa);
	
	Loff = numLBRG.dLoff(aa);
	Boff = numLBRG.dBoff(aa);
	Roff = numLBRG.dRoff(bb);
	
	GTon = numLBRG.dGTon(aa);
	GToff = numLBRG.dGToff(aa);
	
	
	GTpara = [GTon GToff];
	LBRpara = [Lon Bon Ron Loff Boff Roff];
	[stepN OKGO] = RobustGR(Lon,Bon,Ron,Loff,Boff,Roff,dT,vars,glu,GTon,GToff);

	AOK(nn) = OKGO;
	LBRp(nn,:) = LBRpara;
	GTp(nn,:) = GTpara;
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



%====================================================%
%					DATA PREP
%====================================================%
%{
numLBRG = struct('dLon',{PLon},'dBon',{PBon},'dRon',{PRon},...
	'dLoff',{PLoff},'dBoff',{PBoff},'dRoff',{PRoff},...
	'dGTon',{PGTon},'dGToff',{PGToff});
%--		[Lon	Bon		Ron		Loff	Boff	Roff	GTon	GToff]	--%
doLBRG = [0		0		0		0		0		0		1		1	 ];
didLBR = find(doLBRG);


%----------------------------------------------------%
% GT ROBUST TEST

PMX = [GTp AOK' STP'];

METASTEP = [STP' GTp];

STB8 = PMX(((AOK')==1),1:4);	% Stable P-combos (4params)
STB2 = STB8;
DEC8 = PMX(((AOK')==0),1:4);	% Decay  P-combos (4params)
DEC2 = DEC8;
EXP8 = PMX(((AOK')==2),1:4);	% Expand P-combos (4params)
EXP2 = EXP8;


METAS = METASTEP;
METAS(((AOK')==2),1) = Nsteps+(Nsteps-METAS(((AOK')==2),1));





%----------------------------------------------------%
% RATE PARAM ROBUST TEST

PMX = [LBRp AOK' STP'];

didLBR = find(doLBRG);
METASTEP = [STP' LBRp(:,didLBR)];

STB8 = PMX(((AOK')==1),1:8);	% Stable P-combos (8params)
STB2 = STB8(:, didLBR);			% Stable P-combos (2params only)
DEC8 = PMX(((AOK')==0),1:8);	% Decay  P-combos (8params)
DEC2 = DEC8(:, didLBR);			% Decay  P-combos (2params only)
EXP8 = PMX(((AOK')==2),1:8);	% Expand P-combos (8params)
EXP2 = EXP8(:, didLBR);			% Expand P-combos (2params only)


METAS = METASTEP;
METAS(((AOK')==2),1) = Nsteps+(Nsteps-METAS(((AOK')==2),1));
%----------------------------------------------------%
%}
%----------------------------------------------------%
% RATE PARAM ROBUST TEST

LBRGp = [LBRp GTp];

PMX = [LBRGp AOK' STP'];

didLBR = find(doLBRG);


STB10 = PMX(((AOK')==1),1:10);	% Stable P-combos (10params)
STB = STB10(:, didLBR);			% Stable P-combos (2params only)
DEC10 = PMX(((AOK')==0),1:10);	% Decay  P-combos (10params)
DEC = DEC10(:, didLBR);			% Decay  P-combos (2params only)
EXP10 = PMX(((AOK')==2),1:10);	% Expand P-combos (10params)
EXP = EXP10(:, didLBR);			% Expand P-combos (2params only)

METASTEP = [LBRGp(:,didLBR) STP' ];
METAS = METASTEP;
METAS(((AOK')==2),3) = Nsteps+(Nsteps-METAS(((AOK')==2),3));
%----------------------------------------------------%



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
scatter(STB(:,2),STB(:,1), sqrt(STB10(:,10)), 'ob','fill')
ylabel(ylabf1); xlabel(xlabf1);
hold on
scatter(DEC(:,2),DEC(:,1), sqrt(DEC10(:,10)), 'or','fill')
ylabel(ylabf1); xlabel(xlabf1);
hold on
scatter(EXP(:,2),EXP(:,1), sqrt(Nsteps*4), 'og','fill')
ylabel(ylabf1); xlabel(xlabf1);
title('Stable(blue) - Decay(red) - Expand(green)')
%----------------------------------------------------%



%====================================================%
%			FIG2 - 3D Cluster Stability
%----------------------------------------------------%
fh2 = figure(2); hold on
set(fh2,'OuterPosition',(scsz./[2e-3 2e-3 2 1.5]))
%----------------------------------------------------%
xx = METASTEP(:,2);
yy = METASTEP(:,1);
zz = METASTEP(:,3);

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
xx = METAS(:,2);
yy = METAS(:,1);
zz = METAS(:,3);

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
view(280, 40);
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


METASTABLE = LBRp((AOK==1),1:6);
DECAY = LBRp((AOK==0),1:6);
EXPAND = LBRp((AOK==2),1:6);
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




