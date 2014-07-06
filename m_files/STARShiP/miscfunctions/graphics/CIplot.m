function varargout = CIplot(datas,varargin)

dataINo = datas;

if nargin ==1
	cOLOR   = [.01 .95 .01];
	ptype	= 4;
	sbpos   = [.1 .1 .85 .85];
	itemN   = 1;	
elseif nargin ==2
	cOLOR   = varargin{1};
	ptype	= 4;
	sbpos   = [.1 .1 .85 .85];
	itemN   = 1;	
elseif nargin ==3
	cOLOR   = varargin{1};
	ptype	= varargin{2};
	sbpos   = [.1 .1 .85 .85];
	itemN   = 1;
elseif nargin ==4
	cOLOR   = varargin{1};
	ptype	= varargin{2};
	sbpos   = varargin{3};
	itemN   = 1;
elseif nargin ==5
	cOLOR   = varargin{1};
	ptype	= varargin{2};
	sbpos   = varargin{3};
	itemN   = varargin{4};
end




dataIN = dataINo;

if iscell(dataIN)
	datasz = size(dataIN{1});
	if datasz(1) > datasz(2)
		for nn = 1:numel(dataIN)
			dataRE(:,nn) = dataIN{nn};
		end
	else
		for nn = 1:numel(dataIN)
			dataRE(:,nn) = dataIN{nn}';
		end
	end
end


dataIN = {dataRE};


%==============================================%
Mu = mean(dataIN{itemN},2)';
Sd = std(dataIN{itemN},0,2)';
Se = Sd./sqrt(numel(dataIN{itemN}(1,:)));
y_Mu = Mu;
x_Mu = 1:(size(y_Mu,2));
e_Mu = Se;
xx_Mu = 1:0.1:max(x_Mu);
% yy_Mu = spline(x_Mu,y_Mu,xx_Mu);
% ee_Mu = spline(x_Mu,e_Mu,xx_Mu);
yy_Mu = interp1(x_Mu,y_Mu,xx_Mu,'pchip');
ee_Mu = interp1(x_Mu,e_Mu,xx_Mu,'pchip');
p_Mu = polyfit(x_Mu,Mu,3);
x2_Mu = 1:0.1:max(x_Mu);
y2_Mu = polyval(p_Mu,x2_Mu);
%===========================================================%



switch ptype
    case 1 
		subplot('Position',sbpos),...
        confplot(x_Mu,y_Mu,e_Mu);
		hold on;
		hax = gca;
		ph = plot(x_Mu,y_Mu,'o',xx_Mu,yy_Mu);
    case 2
		subplot('Position',sbpos),...
        ph = errorbar(Mu,Se,'x');
		set(ph,'LineStyle','none','LineWidth',1,'Marker','o',...
					'MarkerSize',5,'MarkerEdgeColor',[.2 .2 .2]);
		hax = gca;
    case 3 
		subplot('Position',sbpos),...
        ph = plot(x_Mu,y_Mu,'o',xx_Mu,yy_Mu);
		hax = gca;
	case 4 
		XT_Mu = xx_Mu';
		YT_Mu = yy_Mu';
		ET_Mu = ee_Mu';
		subplot('Position',sbpos),...
        [hl, hp] = boundedline(XT_Mu,YT_Mu, ET_Mu,...
		'cmap',cOLOR(1,:),'alpha','transparency', 0.4);
		ph = hl;
		po = hp;
		hax = gca;
		
	case 5 
		XT_Mu = xx_Mu';
		YT_Mu = yy_Mu';
		ET_Mu = ee_Mu';
		subplot('Position',sbpos),...
        [hl, hp] = boundedline(XT_Mu,YT_Mu, ET_Mu,...
		'cmap',cOLOR(1,:),'alpha','transparency', 0.3);
		%outlinebounds(hl,hp);
		
		XTn=XT_Mu(round(numel(XT_Mu)/2));
		YTn=XT_Mu(round(numel(YT_Mu)/2));
		PLh = light('Position',[XTn YTn 2],'Style','infinite');
		set(hp,'FaceLighting','flat','AmbientStrength',0.8)
		%material dull %shading interp 
		%lighting gouraud; alpha(.4); 
		%lighting phong; alpha(.4);
		ph = hl;
		po = hp;
		hax = gca;
	case 6 
		subplot('Position',sbpos),...
        ciplot((yy_Mu-ee_Mu),(yy_Mu+ee_Mu),xx_Mu,'b');
		% lighting gouraud; alpha(.4); 
		lighting phong; alpha(.4);
		hold on;
		hax = gca;
		ph = plot(x_Mu,y_Mu,'o',xx_Mu,yy_Mu);
    otherwise
        warning('Unexpected plot type.');
end
%===========================================================%








nargchk(0, 2, nargout);
if nargout >= 1
    varargout{1} = ph;
end
if nargout >= 2
    varargout{2} = hax;
end
if nargout >= 3
    
	%varargout{3} = xx_Mu';
	%varargout{3} = yy_Mu';
	%varargout{3} = ee_Mu';
	varargout{3} = Mu;
end

return

%===========================================================%
%===========================================================%







%------------------------------------------%
% [ph2 hax2] = CIplot(data);
% legend([OUTH;ph2],OUTM{:},'Item');
% [LEGH,OBJH,OUTH,OUTM] = legend;
% hold on
%------------------------------------------%
% DATARATE = 100;
% xt = (get(gca,'XTick'))*DATARATE;
% set(gca,'XTickLabel', sprintf('%.0f|',xt))
% set(legend,'Location','NorthWest');
%------------------------------------------%
MS1 = 5; c1=[.9 .1 .9]; c2=[.1 .9 .1];
set(ph,'LineStyle','-','Color',c1,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c1,'MarkerFaceColor',c1);
hTitle  = title('Title');
hXLabel = xlabel('Step');
hYLabel = ylabel('Units (+/- SEM)');
set(gca,'FontName','Helvetica');
set([hTitle, hXLabel, hYLabel],'FontName','AvantGarde');
set([hXLabel, hYLabel],'FontSize',10);
set( hTitle,'FontSize',12,'FontWeight','bold');
set(gca,'Box','off','TickDir','out','TickLength',[.02 .02], ...
'XMinorTick','on','YMinorTick','on','YGrid','on', ...
'XColor',[.3 .3 .3],'YColor',[.3 .3 .3],'LineWidth',1);
% haxes=axis;
% ylim([0 haxes(4)*1.2 ]);
% xlim([0 (haxes(2)*.9)]);





if nargout >= 1
    varargout{1} = ph;
end
if nargout >= 2
    varargout{2} = hax;
end
if nargout >= 3
    varargout{3} = po;
end

return






%===========================================================%
%			ERROR BAR & CONFIDENCE ENVELOPE PLOTS
%===========================================================%
FIGgcf = gcf;
set(FIGgcf,'Units','pixels');scnsize = get(0,'ScreenSize');
pos1 = [scnsize(3)/3  scnsize(4)/5  scnsize(3)/1.2  scnsize(4)/1.7];
set(FIGgcf,'OuterPosition',pos1)
set(gcf,'Color',[.9,.9,.9])

% use: subplot('Position',sbpos),...

%==============================================%
% confplot
% confplot(X,Y,SEM)
%------------------
% plots jagged line and jagged CI envelope
% SEM is distance from Y at each X
% X,Y,SEM must all have same numel
% Y = mean values; X = time; SEM of Y values
%------------------

%{
figure(FIGgcf);
subplot('Position',sbpos),...
confplot(x_Mu,y_Mu,e_Mu);
hold on
%}

figure(1);
confplot(x_Mu,y_Mu,e_Mu);
hold on
hTitle  = title ('confplot');

%==============================================%




%==============================================%
% errorbar
% errorbar(Y,SEM,'ok')
% errorbar(X,Y,SEM,'ok')
%-----------------------
% plots error bars 
% SEM is distance from Y at each X
% X,Y,SEM must all have same numel
% Y = mean values; X = time; SEM of Y values
%-----------------------

%{
subplot('Position',sbpos),...
FIGerb = errorbar(Mu,Se,'x');
set(FIGerb,'LineStyle','none','LineWidth',1,'Marker','o',...
'MarkerSize',5,'MarkerEdgeColor',[.2 .2 .2]);
hold on
%}

figure(1);
FIGerb = errorbar(Mu,Se,'x');
set(FIGerb,'LineStyle','none','LineWidth',1,'Marker','o',...
'MarkerSize',5,'MarkerEdgeColor',[.2 .2 .2]);
hTitle  = title ('errorbar');



%==============================================%




%==============================================%
% spline
%-----------------------
% Smoothed Plot Line (SPLine)
%
% 1) create fine grain x-axis timesteps
%	 >> X=timesteps; Y=Mu;
%    >> XX = 1:0.1:max(X);
%
% 2) use spline(X,Y,XX)
%    >> YY = spline(X,Y,XX);
%
% 3) plot smooth line
%	 >> plot(X,Y,'o',XX,YY);
% 
%-----------------------

%{
subplot('Position',sbpos),...
plot(x_Mu,y_Mu,'o',xx_Mu,yy_Mu);
hold on
%}

figure(1);
plot(x_Mu,y_Mu,'o',xx_Mu,yy_Mu);
hTitle  = title ('spline');

%==============================================%



%==============================================%
% boundedline
% boundedline(xx, yy, ee)
%-----------------------
% plots boundedline
%-----------------------

%{
subplot('Position',sbpos),...
[hl, hp] = boundedline(xx_Mu, yy_Mu, ee_Mu);
hold on
lighting gouraud; 
alpha(.4)


% XT_Mu = xx_Mu';
% YT_Mu = yy_Mu';
% ET_Mu = ee_Mu';
% %XT_Mu(2,:) = sqrt(xx_Mu)
% YT_Mu(2,:) = YT_Mu-(YT_Mu./2)
% ET_Mu(2,:) = ET_Mu-(ET_Mu./2)
% XT_Mu = XT_Mu';
% YT_Mu = YT_Mu';
% ET_Mu = ET_Mu';
% figure(2);
% [hl, hp] = boundedline(XT_Mu,YT_Mu, ET_Mu,...
% 'cmap', [.4 .1 .3; .1 .8 .7],'alpha','transparency', 0.07);
% 'cmap',lines(4),'alpha','transparency', 0.07);
% lighting gouraud; alpha(.4);
%}

figure(1);
[hl, hp] = boundedline(xx_Mu, yy_Mu, ee_Mu);
hTitle  = title ('boundedline');
lighting gouraud; 
alpha(.4)

%==============================================%


%==============================================%
% ciplot
% ciplot((YY-EE),(YY+EE),XX,'b')
%-----------------------
% confidence interval plot (CIPlot)
%
% 1) create fine grain x-axis timesteps
%	 >> X=timesteps; Y=Mu;
%    >> XX = 1:0.1:max(X);
%
% 2) use spline function to get YY and EE
%    >> YY = spline(X,Y,XX);
%    >> EE = spline(X,Y,XX);
%
% 3) plot
%	 >> ciplot((YY-EE),(YY+EE),XX,'b')
% 
%-----------------------

%{
subplot('Position',sbpos),...
ciplot((yy_Mu-ee_Mu),(yy_Mu+ee_Mu),xx_Mu,'b')
%ciplot((y_Mu-e_Mu),(y_Mu+e_Mu),x_Mu,'b')
% lighting gouraud; alpha(.4);
%}

figure(1);
ciplot((yy_Mu-ee_Mu),(yy_Mu+ee_Mu),xx_Mu,'b')
hTitle  = title ('ciplot');
lighting gouraud; 
alpha(.4)



%==============================================%






varargout={dataIN};


end








