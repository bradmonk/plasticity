clc, clear all, close all


Boff = 0;

for n=1:10
	
Nsteps = 80000;
dT = .01;
Lon = 1.5;
Loff = 1.5;
Bon = 0;
% Boff = 1;
Ron = .01;
Roff = 1;

Boff = Boff+1;

QRvals = [Bon Ron Boff Roff Lon Loff];
params = [Nsteps dT];


	
for m=1:10
[NSteps,stepN,S0,S] = ShouvalORMFUN(QRvals,params);

% NStepsOUT(m) = Nsteps;
stepNOUT(m) = stepN;
% S0OUT{m} = S0;
% SOUT{m} = S;
% CSenOUT{m} = CSen;
% CSexOUT{m} = CSex;
% NoffNOUT(m) = NoffN;
% SoffsOUT{m} = Soffs;
% SonsOUT{m} = Sons;
% doNoffOUT(m) = doNoff;

end

SLife(n,:) = stepNOUT;

end

SSLife = SLife.*dT;
SLF = {log10(SSLife)};


% SHOUVAL'S WAY
%{.
Boff = 0;
for m = 1:10
Boff=Boff+1;
NBRS = 1;
Pkoff = (2 .* Roff .* exp(-Boff * (NBRS-Loff))) ./ (1 + exp(-Boff * (NBRS-Loff))) .* dT;
LF1Time(m) = dT / Pkoff;
end

Boff = 0;
for m = 1:10
Boff=Boff+1;
NBRS = 2;
Pkoff = (2 .* Roff .* exp(-Boff * (NBRS-Loff))) ./ (1 + exp(-Boff * (NBRS-Loff))) .* dT;
LF2Time(m) = dT / Pkoff;
end

Boff = 0;
for m = 1:10
Boff=Boff+1;
NBRS = 3;
Pkoff = (2 .* Roff .* exp(-Boff * (NBRS-Loff))) ./ (1 + exp(-Boff * (NBRS-Loff))) .* dT;
LF3Time(m) = dT / Pkoff;
end

Boff = 0;
for m = 1:10
Boff=Boff+1;
NBRS = 4;
Pkoff = (2 .* Roff .* exp(-Boff * (NBRS-Loff))) ./ (1 + exp(-Boff * (NBRS-Loff))) .* dT;
LF4Time(m) = dT / Pkoff;
end

%}

% BRADS WAY
%{
Boff = 0;
for m = 1:10
Boff=Boff+1;
NBRS = 1;
Lhoff = ((-NBRS)+Loff) .* (-Boff);
Poff = 1 ./ (1+exp(Lhoff));
Pkoff = ( Roff * dT * Poff );
LF1Time(m) = dT / Pkoff;
end

Boff = 0;
for m = 1:10
Boff=Boff+1;
NBRS = 2;
Lhoff = ((-NBRS)+Loff) .* (-Boff);
Poff = 1 ./ (1+exp(Lhoff));
Pkoff = ( Roff * dT * Poff );
LF2Time(m) = dT / Pkoff;
end

Boff = 0;
for m = 1:10
Boff=Boff+1;
NBRS = 3;
Lhoff = ((-NBRS)+Loff) .* (-Boff);
Poff = 1 ./ (1+exp(Lhoff));
Pkoff = ( Roff * dT * Poff );
LF3Time(m) = dT / Pkoff;
end

Boff = 0;
for m = 1:10
Boff=Boff+1;
NBRS = 4;
Lhoff = ((-NBRS)+Loff) .* (-Boff);
Poff = 1 ./ (1+exp(Lhoff));
Pkoff = ( Roff * dT * Poff );
LF4Time(m) = dT / Pkoff;
end
%}

LogLF1Time = [];
LogLF2Time = [];
LogLF3Time = [];
LogLF4Time = [];

LogLF1Time = abs(log10(LF1Time));
LogLF2Time = log10(LF2Time);
LogLF3Time = log10(LF3Time);
LogLF4Time = log10(LF4Time);

LoogLF1Time = [LogLF1Time;(LogLF1Time.*.99);(LogLF1Time.*1.01);LogLF1Time;(LogLF1Time.*.99)];
LoogLF2Time = [LogLF2Time;(LogLF2Time.*.99);(LogLF2Time.*1.01);LogLF2Time;(LogLF2Time.*.99)];
LoogLF3Time = [LogLF3Time;(LogLF3Time.*.99);(LogLF3Time.*1.01);LogLF3Time;(LogLF3Time.*.99)];
LoogLF4Time = [LogLF4Time;(LogLF4Time.*.99);(LogLF4Time.*1.01);LogLF4Time;LogLF4Time];

LgLF1Time = {LoogLF1Time'};
LgLF2Time = {LoogLF2Time'};
LgLF3Time = {LoogLF3Time'};
LgLF4Time = {LoogLF4Time'};

%======================================================================%
%							FIGURE SETUP
%----------------------------------------------------------------------%
fig21 = figure(21);
set(21,'Units','pixels');scnsize = get(0,'ScreenSize');
pos1 = [scnsize(3)/3  scnsize(4)/5  scnsize(3)/3.1  scnsize(4)/1.5];
set(fig21,'OuterPosition',pos1)
set(gcf,'Color',[.9,.9,.9])
%----------------------------------------------------------------------%
c1= [.9 .2 .2]; c2= [.2 .4 .6]; c3= [.4 .8 .4]; c4= [.6 .6 .6];
c11=[.9 .3 .3]; c22=[.3 .5 .7]; c33=[.5 .9 .5]; c44=[.7 .7 .7];
%===========================================================%
% FIG1 TOP LEFT: GluR Subtypes in Synapses
%===========================================================%
cOLOR = [c1; c2; c3; c4; c11; c22; c33; c44];
sbpos = [.1 .1 .8 .8]; ptype = 4;
%---
itemN = 1;  
[ph1 hax1] = CIplot(SLF,sbpos,itemN,cOLOR,ptype);
leg1 = legend(ph1,'Cluster Sim');
[LEGH,OBJH,OUTH,OUTM] = legend;
hold on
%---
itemN = 1;  
[ph2 hax2] = CIplot(LgLF1Time,sbpos,itemN,cOLOR,ptype);
legend([OUTH;ph2],OUTM{:},'1-Neighbor');
[LEGH,OBJH,OUTH,OUTM] = legend;
hold on
%---
itemN = 1;
[ph3 hax3] = CIplot(LgLF2Time,sbpos,itemN,cOLOR,ptype);
legend([OUTH;ph3],OUTM{:},'2-Neighbors');
[LEGH,OBJH,OUTH,OUTM] = legend;
hold on
%---
itemN = 1;
[ph4 hax4] = CIplot(LgLF3Time,sbpos,itemN,cOLOR,ptype);
legend([OUTH;ph4],OUTM{:},'3-Neighbors');
[LEGH,OBJH,OUTH,OUTM] = legend;
hold on
%---
itemN = 1;
[ph5 hax5] = CIplot(LgLF4Time,sbpos,itemN,cOLOR,ptype);
legend([OUTH;ph5],OUTM{:},'4-Neighbors');
[LEGH,OBJH,OUTH,OUTM] = legend;
hold on
%---
%------------------------------------------%
set(legend,'Location','NorthWest');
set(legend,'FontSize',16);
legend('off')
%------------------------------------------%
MS1 = 7; MS2 = 2; clr1=[.1 .1 .1];clr2=[.1 .1 .1];
set(ph1,'LineStyle','-','Color',clr1,'LineWidth',1,...
'Marker','s','MarkerSize',7,'MarkerEdgeColor',clr1,'MarkerFaceColor',clr1);
set(ph2,'LineStyle',':','Color',c2,'LineWidth',4,...
'Marker','none','MarkerSize',MS1,'MarkerEdgeColor',c2,'MarkerFaceColor',c22);
set(ph3,'LineStyle','-.','Color',c3,'LineWidth',4,...
'Marker','none','MarkerSize',MS1,'MarkerEdgeColor',c3,'MarkerFaceColor',c33);
set(ph4,'LineStyle','-','Color',c4,'LineWidth',4,...
'Marker','none','MarkerSize',MS1,'MarkerEdgeColor',c4,'MarkerFaceColor',c44);
set(ph5,'LineStyle','--','Color',c1,'LineWidth',3,...
'Marker','none','MarkerSize',MS1,'MarkerEdgeColor',c1,'MarkerFaceColor',c11);
hTitle  = title('Cluster vs. Particle Lifetimes');
hXLabel = xlabel('Beta');
hYLabel = ylabel('Lifetime (log steps +/- SEM) ');
set(gca,'FontName','Century Gothic');
set([hTitle, hXLabel, hYLabel],'FontName','Arial');
set([hXLabel, hYLabel],'FontSize',16);
% set( hTitle,'FontSize',18,'FontWeight','bold');
set(hTitle,'FontSize',18);
set(gca,'Box','off','TickDir','out','TickLength',[.02 .02], ...
'XMinorTick','on','YMinorTick','on','YGrid','on', ...
'XColor',[.3 .3 .3],'YColor',[.3 .3 .3],'LineWidth',1);
% haxes=axis;
% ylim([0 haxes(4)*1.2 ]);
% xlim([0 (haxes(2))]);
%======================================================================%




%======================================================================%
% FONT NAME TESTERS
%{
%----------------------------------------------------------------------%
fig21 = figure(21);
set(21,'Units','pixels');scnsize = get(0,'ScreenSize');
pos1 = [scnsize(3)/3  scnsize(4)/5  scnsize(3)/1.5  scnsize(4)/1.5];
set(fig21,'OuterPosition',pos1)
set(gcf,'Color',[.9,.9,.9])
%------


hTXT = text(.1,.9,0,'Century Gothic Font Name');
set(hTXT,'FontName','Century Gothic','FontSize',18);
% set(hTXT,'FontWeight','bold');

hTXT = text(.1,.8,0,'Book Antiqua Font Name');
set(hTXT,'FontName','Book Antiqua','FontSize',18);
% set(hTXT,'FontWeight','bold');

hTXT = text(.1,.7,0,'Arial Font Name');
set(hTXT,'FontName','Arial','FontSize',18);
% set(hTXT,'FontWeight','bold');

hTXT = text(.1,.6,0,'Microsoft Sans Serif Font Name');
set(hTXT,'FontName','Microsoft Sans Serif','FontSize',18);
% set(hTXT,'FontWeight','bold');

hTXT = text(.1,.5,0,'Calibri Font Name');
set(hTXT,'FontName','Calibri','FontSize',18);
% set(hTXT,'FontWeight','bold');

hTXT = text(.1,.4,0,'Times New Roman Font Name');
set(hTXT,'FontName','Times New Roman','FontSize',18);
% set(hTXT,'FontWeight','bold');

hTXT = text(.1,.3,0,'Verdana Font Name');
set(hTXT,'FontName','Verdana','FontSize',18);
% set(hTXT,'FontWeight','bold');

hTXT = text(.1,.2,0,'xxxxxx Font Name');
set(hTXT,'FontName','xxxxxx','FontSize',18);
% set(hTXT,'FontWeight','bold');

hTXT = text(.1,.1,0,'xxxxxx TT Font Name');
set(hTXT,'FontName','xxxxxx','FontSize',18);
% set(hTXT,'FontWeight','bold');
%}
%----------------------------------------------------------------------%


