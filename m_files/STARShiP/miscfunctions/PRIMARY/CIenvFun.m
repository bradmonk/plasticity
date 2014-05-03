function varargout = CIenvFun(varargin)


% if nargin('CIenvFun') < 1
% 	sbpos = [.05 .09 .43 .8];
% 	itemN = 1; 
% 	c1= [.9 .2 .2]; c2= [.2 .4 .6]; c3= [.4 .8 .4]; c4= [.6 .6 .6];
% 	c11=[.9 .3 .3]; c22=[.3 .5 .7]; c33=[.5 .9 .5]; c44=[.7 .7 .7];
% 	cOLOR = [c1; c2; c3; c4];
% 	
% 	for m = 1:4
% 		for n = 1:4
% 			pholder2(:,n) = abs(randn(10,1)) .* ((1:10).^2)';
% 		end
% 		pholder3{m}=pholder2;
% 	end
% 
% varargset = {pholder3,sbpos,itemN,cOLOR};
% 
% 
% dataINo = varargset{1};
% sbpos   = varargset{2};
% itemN   = varargset{3};
% cOLOR   = varargset{4};
% else
% %===========================================================%
% % varargin = {pholder3,sbpos,itemN};
% dataINo = varargin{1};
% sbpos   = varargin{2};
% itemN   = varargin{3};
% cOLOR   = varargin{4};
% end

dataINo = varargin{1};
sbpos   = varargin{2};
itemN   = varargin{3};
cOLOR   = varargin{4};
ptype	= varargin{5};

dataTemp = [];
dataIN = [];

% for c = 1:numel(dataINo); for now = 1:numel(dataINo{1}(1,:));
% dataTemp(:,now) = dataINo{1,now}(:,c);
% end; dataIN{c}=dataTemp; end;


dataIN = [];
dataTemp = [];
for c = 1:numel(dataINo{1}(1,:));
	for now = 1:numel(dataINo); 
		dataTemp(:,now) = dataINo{1,now}(:,c);
	end; 
dataIN{c}=dataTemp; 
end;


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



%{

% for c = 1:numel(redater); for now = 1:numel(redater{itemN}(1,:));
% redater2(:,now) = redater{1,now}(:,c);
% end; redater3{c}=redater2; end;


% reda3Mu = [mean(redater3{1},2)]';
% reda3Sd = [std(redater3{1},0,2)]';
% reda3Se = [reda3Sd./sqrt(numel(redater3{1}(1,:)))];
% 
% t_reda3Mu = 1:(size(reda3Mu,2));
% p_reda3Mu = polyfit(t_reda3Mu,reda3Mu,6);
% t2_reda3Mu = 1:0.1:max(t_reda3Mu);
% y2_reda3Mu = polyval(p_reda3Mu,t2_reda3Mu);
% %[y2p,y2d] = polyconf(p_reda3Mu,t_reda3Mu,y2_reda3Mu)
% 
% y_Mu = reda3Mu;
% x_Mu = 1:(size(y_Mu,2));
% e_Mu = reda3Se;
% xx_Mu = 1:0.1:max(x_Mu);
% yy_Mu = spline(x_Mu,y_Mu,xx_Mu);


% resSe1 = [];
% linstps = round(numel(yy_Mu)/numel(reda3Se));
% redo=1;
% for red2 = 2:numel(reda3Se)
% red1=red2-1;
% resSe1(:,redo:(redo-1+linstps)) = linspace(reda3Se(red1),reda3Se(red2),linstps);
% redo = redo+linstps;
% end
% resSe1(numel(resSe1)+1) = resSe1(numel(resSe1));
% 
% ee_Mu = resSe1;
%}
%{
% plottype = 'confplot';
% plottype = 'errorbar';
% plottype = 'spline';
% plottype = 'boundedline';
% plottype = 'ciplot';

switch plottype
    case 'confplot' 
        confplot(x_Mu,y_Mu,e_Mu);
    case 'errorbar'
        ebf = errorbar(Mu,Se,'x');
		set(ebf,'LineStyle','none','LineWidth',1,'Marker','o',...
					'MarkerSize',5,'MarkerEdgeColor',[.2 .2 .2]);
    case 'spline' 
        plot(x_Mu,y_Mu,'o',xx_Mu,yy_Mu);
	case 'boundedline' 
        [hl, hp] = boundedline(xx_Mu, yy_Mu, ee_Mu);
		lighting gouraud; alpha(.4);
	case 'ciplot' 
        ciplot((yy_Mu-ee_Mu),(yy_Mu+ee_Mu),xx_Mu,'b');
		lighting gouraud; alpha(.4);
    otherwise
        warning('Unexpected plot type.');
end
%}
%{
%===========================================================%
% 1 = 'confplot';
% 2 = 'errorbar';
% 3 = 'spline';
% 4 = 'boundedline';
% 5 = 'ciplot';
plottype = 4;

switch plottype
    case 1 
        confplot(x_Mu,y_Mu,e_Mu);
    case 2
        ebf = errorbar(Mu,Se,'x');
		set(ebf,'LineStyle','none','LineWidth',1,'Marker','o',...
					'MarkerSize',5,'MarkerEdgeColor',[.2 .2 .2]);
    case 3 
        plot(x_Mu,y_Mu,'o',xx_Mu,yy_Mu);
	case 4 
        [hl, hp] = boundedline(xx_Mu, yy_Mu, ee_Mu,...
			'alpha','cmap', cool(4), 'transparency', 0.5);
		lighting gouraud; alpha(.4);
	case 5 
        ciplot((yy_Mu-ee_Mu),(yy_Mu+ee_Mu),xx_Mu,'b');
		lighting gouraud; alpha(.4);
    otherwise
        warning('Unexpected plot type.');
end
%===========================================================%
%}


%{.
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
%}


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








