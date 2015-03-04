%----------------------------------------------------------------------
% Modified by Amie Adams
%Adapted by Bruce Land from:
% Physiology 317 - Methods in Computational Neurobiology
% HODGKIN-HUXLEY: Response of Squid Axon to Current Injection
% M. Nelson  1995, 1997
%----------------------------------------------------------------------

clear
clf
orient landscape

%%%%%%%%%%%%%%%%%%%%%%%
% INITIALIZATION
%%%%%%%%%%%%%%%%%%%%%%%
DT = 0.025;
TMAX = 50;
t = 0:DT:TMAX;

VREST = -60.0;

GNA = 120.0;        % mS/cm^2
GK = 36.0;        % mS/cm^2
GLEAK = 1.0;        % mS/cm^2 0.3

ENA = 115.0 + VREST;        % mV
EK = -12.0 + VREST;        % mV
ELEAK = 10.613 + VREST;        % mV

C_M = 1.0;        % uF/cm^2

%set the axon properties
specR = 30 ; % ohm cm
dia = 0.005 ; % cm 0.05
nodeDist = 0.005 ; %cm
Raxon = specR*nodeDist/(pi*(dia^2)/4) ; % R=specR*length/area
gAxon = 1000/Raxon ; %invert and convert to milliS
nNodes = 20 ;
nTime = length(t) ;

%current strength
INJECT_LIST = 30;


%%%%%%%%%%%%%%%%%%%%%%%
% MAIN SIMULATION LOOP
%%%%%%%%%%%%%%%%%%%%%%%

clear v  m  h  n                                % clear old varibles
v = ones(nTime, nNodes)*VREST; % initial membrane voltage
m = ones(nTime,nNodes).*alpha_m(v-VREST)./(alpha_m(v-VREST)+beta_m(v-VREST));        %initial (steady-state) m
h = ones(nTime,nNodes).*alpha_h(v-VREST)./(alpha_h(v-VREST)+beta_h(v-VREST));        %initial (steady-state) h
n = ones(nTime,nNodes).*alpha_n(v-VREST)./(alpha_n(v-VREST)+beta_n(v-VREST));        %initial (steady-state) n

for i=2:length(t)
    I = [INJECT_LIST*(t(i)>1 & t(i)<=2), zeros(1,nNodes-1)] ;
    M = m(i-1,:);        % get values from last time step
    H = h(i-1,:);        % (hopefully this substitution makes
    N = n(i-1,:);        % the following code a bit easier to read)
    V = v(i-1,:);

    gNa = GNA * M.^3 .* H ;
    gK  = GK  * N.^4;

    mdot = alpha_m(V-VREST).*(1-M) - beta_m(V-VREST).*M;
    hdot = alpha_h(V-VREST).*(1-H) - beta_h(V-VREST).*H;
    ndot = alpha_n(V-VREST).*(1-N) - beta_n(V-VREST).*N;
    vdot = (I - gNa.*(V-ENA) - gK.*(V-EK) - GLEAK.*(V-ELEAK))/C_M;

    for j=2:nNodes-1
        vdot(j) = vdot(j) + gAxon * (V(j-1)-V(j)-V(j)+V(j+1)) ;
    end
    vdot(1) = vdot(1) - gAxon * (V(2)-V(1)) ;
    vdot(nNodes) = vdot(nNodes) + gAxon * (V(nNodes-1)-V(nNodes)) ;

    m(i,:) = m(i-1,:) + mdot*DT;                % Euler integration
    h(i,:) = h(i-1,:) + hdot*DT;
    n(i,:) = n(i-1,:) + ndot*DT;
    v(i,:) = v(i-1,:) + vdot*DT;
end

figure(1);clf;
cmap = jet(nNodes);
%plot(t,v)
for j=1:nNodes
   % subplot(nNodes,1,j);
    plot(t,v(:,j),'color',cmap(j,:));
    hold on
    ylabel('v (mV)')
    %xlabel('t (msec)'); ylabel('v (mV)')
    set(gca,'ylim',[-100 50])
end

%get speed by finding peaks for each trace
for j=1:nNodes
    [vPeak(j),index(j)] = max(v(:,j));
end
%m/sec = (cm/100)/(mSec/1000)
10*nodeDist./(t(index(end-2))-t(index(end-3)))

figure(2);clf;
%plot peak arrival time versus distance
x = [0:nNodes-1]*nodeDist;
y = t(index);
plot(x, y,'x')
hold on

%least squares linear fit
p = polyfit(x,y,1);
plot(x,polyval(p,x),'r')
distance=(x(end)-x(1));
time=(y(end)-y(1));
%print velocity
title(['velocity=',num2str(10*distance/time,'%4.4f'),' meters/sec'])

