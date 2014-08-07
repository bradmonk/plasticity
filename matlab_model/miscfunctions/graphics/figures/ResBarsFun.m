function varargout = ResBarsFun(varargin)




%%
%==========================================================%
% For example, suppose you measure a quantity y at several values of time t:

t = [0 0.3 0.8 1.1 1.6 2.3];
y = [0.6 0.67 1.01 1.35 1.47 1.25];
plot(t,y,'o')


% You can try modeling this data using a second-degree polynomial function:
% The unknown coefficients a0, a1, and a2 are computed by minimizing the sum of the 
% squares of the deviations of the data from the model (least-squares fit).
% To find the polynomial coefficients, type the following at the MATLAB prompt:

p=polyfit(t,y,2)


% MATLAB calculates the polynomial coefficients in descending powers:
% p = -0.2942    1.0231    0.4981
% The second-degree polynomial model of the data is given by the following equation:
% To plot the model with the data, evaluate the polynomial at uniformly spaced times 
% t2 and overlay the original data on a plot:

t2 = 0:0.1:2.8;     % Define a uniformly spaced time vector
y2=polyval(p,t2);   % Evaluate the polynomial at t2
figure
plot(t,y,'o',t2,y2) % Plot the fit on top of the data 
                    % in a new Figure window


					
% Use the following syntax to calculate the residuals:

y2=polyval(p,t); % Evaluate model at the data time vector
res=y-y2; % Calculate the residuals by subtracting
figure, plot(t,res,'+') % Plot the residuals


% Repeat the exercise, this time using a fifth-degree polynomial from polyfit:

p5= polyfit(t,y,5)

% p5 = 0.7303   -3.5892    5.4281   -2.5175    0.5910    0.6000

y3 = polyval(p5,t2);   % Evaluate the polynomial at t2
figure
plot(t,y,'o',t2,y3) % Plot the fit on top of the data 
                    % in a new Figure window
%%
%==========================================================%

% errorbar(Y,E) plots Y and draws an error bar at each element of Y. The error bar 
% is a distance of E(i) above and below the curve so that each bar is symmetric and 2*E(i) long.

%==========================================================%
MU1 = [1 2];
SIGMA1 = [2 0; 0 .5];
MU2 = [-3 -5];
SIGMA2 = [1 0; 0 1];
X = [mvnrnd(MU1,SIGMA1,1000);
mvnrnd(MU2,SIGMA2,1000)];
scatter(X(:,1),X(:,2),10,'.')
options = statset('Display','final');
obj = gmdistribution.fit(X,2,'Options',options);
hold on
h = ezcontour(@(x,y)pdf(obj,[x y]),[-8 6],[-8 6]);
hold off
%==========================================================%

y1 = [2 3 5 7 8 5];

t = 1:(size(y,2));
p=polyfit(t,y,2)
t2 = 0:0.1:max(y);     % Define a uniformly spaced time vector
y2=polyval(p,t2);   % Evaluate the polynomial at t2
figure
plot(t,y,'o',t2,y2)
errorbar(y,E)



%%
%==========================================================%
redaterO=reDATAG2SLOTSDATA;
redater = redaterO;

for c = 1:numel(redater); for now = 1:numel(redater{1}(1,:));
redater2(:,now) = redater{1,now}(:,c);
end; redater3{c}=redater2; end;

reda3Mu = [mean(redater3{1},2)]';
reda3Sd = [std(redater3{1},0,2)]';
reda3Se = [reda3Sd./sqrt(numel(redater3{1}(1,:)))];

t_reda3Mu = 1:(size(reda3Mu,2));
p_reda3Mu = polyfit(t_reda3Mu,reda3Mu,6);
t2_reda3Mu = 1:0.1:max(t_reda3Mu);
y2_reda3Mu = polyval(p_reda3Mu,t2_reda3Mu);
%[y2p,y2d] = polyconf(p_reda3Mu,t_reda3Mu,y2_reda3Mu)

y_Mu = reda3Mu;
x_Mu = 1:(size(y_Mu,2));
xx_Mu = 1:0.1:max(x_Mu);
yy_Mu = spline(x_Mu,y_Mu,xx_Mu);


resSe1 = [];
linstps = round(numel(yy_Mu)/numel(reda3Se));
redo=1;
for red2 = 2:numel(reda3Se)
red1=red2-1;
resSe1(:,redo:(redo-1+linstps)) = linspace(reda3Se(red1),reda3Se(red2),linstps);
redo = redo+linstps;
end
resSe1(numel(resSe1)+1) = resSe1(numel(resSe1));

gMu1 = figure(1);
% geMu1 = errorbar(reda3Mu,reda3Se,'x');
% set(geMu1,'LineStyle','none','LineWidth',1,'Marker','o',...
% 'MarkerSize',5,'MarkerEdgeColor',[.2 .2 .2]);
hold on
gpMu1 = plot(x_Mu,y_Mu,'o',xx_Mu,yy_Mu);
set(gpMu1,'Color',[0 0 .5],'LineWidth',1,'Marker','.');
[hl, hp] = boundedline(xx_Mu, yy_Mu, resSe1);





%%
%==========================================================%


%%
%==========================================================%
redaterO=reDATAG2SLOTSDATA;
redater = redaterO;

for c = 1:numel(redater); for now = 1:numel(redater{1}(1,:));
redater2(:,now) = redater{1,now}(:,c);
end; redater3{c}=redater2; end;

reda3Mu = [mean(redater3{1},2)]';
reda3Sd = [std(redater3{1},0,2)]';
reda3Se = [reda3Sd./sqrt(numel(redater3{1}(1,:)))];

t_reda3Mu = 1:(size(reda3Mu,2));
p_reda3Mu = polyfit(t_reda3Mu,reda3Mu,6);
t2_reda3Mu = 1:0.1:max(t_reda3Mu);
y2_reda3Mu = polyval(p_reda3Mu,t2_reda3Mu);
%[y2p,y2d] = polyconf(p_reda3Mu,t_reda3Mu,y2_reda3Mu)

y_Mu = reda3Mu;
x_Mu = 1:(size(y_Mu,2));
e_Mu = reda3Se;
xx_Mu = 1:0.1:max(x_Mu);
yy_Mu = spline(x_Mu,y_Mu,xx_Mu);


resSe1 = [];
linstps = round(numel(yy_Mu)/numel(reda3Se));
redo=1;
for red2 = 2:numel(reda3Se)
red1=red2-1;
resSe1(:,redo:(redo-1+linstps)) = linspace(reda3Se(red1),reda3Se(red2),linstps);
redo = redo+linstps;
end
resSe1(numel(resSe1)+1) = resSe1(numel(resSe1));

gMu1 = figure(1);
geMu1 = errorbar(reda3Mu,reda3Se,'x');
set(geMu1,'LineStyle','none','LineWidth',1,'Marker','o',...
'MarkerSize',5,'MarkerEdgeColor',[.2 .2 .2]);
confplot(x_Mu,y_Mu,e_Mu)
hold on
% gpMu1 = plot(x_Mu,y_Mu,'o',xx_Mu,yy_Mu);
% set(gpMu1,'Color',[0 0 .5],'LineWidth',1,'Marker','.');
% [hl, hp] = boundedline(xx_Mu, yy_Mu, resSe1);





%%
%==========================================================%

mu = @(x) -1.9+.23*x;
x = 5:.1:15;
yhat = mu(x);
dy = -3.5:.1:3.5; sz = size(dy); k = (length(dy)+1)/2;
x1 =  7*ones(sz); y1 = mu(x1)+dy; z1 = normpdf(y1,mu(x1),1);
x2 = 10*ones(sz); y2 = mu(x2)+dy; z2 = normpdf(y2,mu(x2),1);
x3 = 13*ones(sz); y3 = mu(x3)+dy; z3 = normpdf(y3,mu(x3),1);
plot3(x,yhat,zeros(size(x)),'b-', ...
      x1,y1,z1,'r-', x1([k k]),y1([k k]),[0 z1(k)],'r:', ...
      x2,y2,z2,'r-', x2([k k]),y2([k k]),[0 z2(k)],'r:', ...
      x3,y3,z3,'r-', x3([k k]),y3([k k]),[0 z3(k)],'r:');
zlim([0 1]);
xlabel('X'); ylabel('Y'); zlabel('Probability density');
grid on; view([-45 45]);

%%
%==========================================================%

mu = @(x) exp(-1.9+.23*x);
x = 5:.1:15;
yhat = mu(x);
x1 =  7*ones(1,5);  y1 = 0:4; z1 = poisspdf(y1,mu(x1));
x2 = 10*ones(1,7); y2 = 0:6; z2 = poisspdf(y2,mu(x2));
x3 = 13*ones(1,9); y3 = 0:8; z3 = poisspdf(y3,mu(x3));
plot3(x,yhat,zeros(size(x)),'b-', ...
      [x1; x1],[y1; y1],[z1; zeros(size(y1))],'r-', x1,y1,z1,'r.', ...
      [x2; x2],[y2; y2],[z2; zeros(size(y2))],'r-', x2,y2,z2,'r.', ...
      [x3; x3],[y3; y3],[z3; zeros(size(y3))],'r-', x3,y3,z3,'r.');
zlim([0 1]);
xlabel('X'); ylabel('Y'); zlabel('Probability');
grid on; view([-45 45]);





%%
%==========================================================%








%%
%==========================================================%


end