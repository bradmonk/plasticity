function varargout = CIplot(varargin)


dataINo = varargin{1};
sbpos   = varargin{2};
itemN   = varargin{3};
cOLOR   = varargin{4};
ptype	= varargin{5};



dataIN = dataINo;

%==============================================%
Mu = [mean(dataIN{itemN},2)]';
Sd = [std(dataIN{itemN},0,2)]';
Se = [Sd./sqrt(numel(dataIN{itemN}(1,:)))];
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
plottype = ptype;

switch plottype
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
		'cmap',cOLOR(itemN,:),'alpha','transparency', 0.4);
		ph = hl;
		po = hp;
		hax = gca;
	case 5 
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
    varargout{3} = po;
end
% varargout={dataIN;hl};
% varargout={hl;hp};
return
%===========================================================%
%===========================================================%




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








