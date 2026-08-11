function [Nsteps,S1pop2_num,S2beta,S2tau,SC2beta,SCr,initial_up,sz1_bound,...
S,S1r,S2delta_t,SC1L,SC2deltaT,SCro,initial_up1,sz2,...
S1,S1ro,S2mem_num,SC1L2,SC2mu,SCszb,initial_up2,sz2_bound,...
S1L1,S1ro2,S2mu,SC1beta,SC2r,SCszi,sap1,sz_bound,...
S1L1_array,S1s,S2mu_array,SC1deltaT,SC2ro,SCtau,sap2,time,...
S1L2,S1sum,S2num_epochs,SC1mu,SC2szb,center,sap3,...
S1L2_array,S1sumOrig,S2pop2_num,SC1r,SC2szi,center1,sap4,...
S1beta,S1tau,S2r,SC1ro,SC2tau,center2,sap5,...
S1delta_t,S2,S2ro,SC1szb,SCL,h_mask,sap6,...
S1mem_num,S2L1,S2ro2,SC1szi,SCL2,initial_down,sap7,...
S1mu,S2L1_array,S2s,SC1tau,SCbeta,initial_down1,sap8,...
S1mu_array,S2L2,S2sum,SC2L,SCdeltaT,initial_down2,sz,...
S1num_epochs,S2L2_array,S2sumOrig,SC2L2,SCmu,initial_sz,sz1] = SClusters(Nsteps,...
sap1, sap2, sap3, sap6, sap4, sap7, sap5, sap8,runSAPPSD1,runSAPPSD2);

%{
SCdeltaT = 0.01;	% Shouval [.01]		Brad [.01]
SCbeta = 60;		% Shouval [60]		Brad [50]
SCtau = 1.0;		% Shouval [1.0]		Brad [1.8]
SCmu = 1/SCtau;		% Shouval [1/tau]	Brad [1/tau]
SCL = 1.5;			% Shouval [1.5]		Brad [1.2]
SCL2 = 1.1;			% Shouval [0.9]		Brad [1.1]
SCr = 10;			% Shouval [10]		Brad [15]
SCro = 0.95;		% Shouval [.95]		Brad [.90]
SCszi = 7;			% Shouval [8]		Brad [7]
SCszb = 17;			% Shouval [17]		Brad [17]

SCdeltaT = 0.01;	% Shouval [.01]		Brad [.01]
SCbeta = 50;		% Shouval [60]		Brad [50]
SCtau = 1.8;		% Shouval [1.0]		Brad [1.8]
SCmu = 1/SCtau;		% Shouval [1/tau]	Brad [1/tau]
SCL = 1.2;			% Shouval [1.5]		Brad [1.2]
SCL2 = 1.1;			% Shouval [0.9]		Brad [1.1]
SCr = 15;			% Shouval [10]		Brad [15]
SCro = 0.90;		% Shouval [.95]		Brad [.90]
SCszi = 7;			% Shouval [8]		Brad [7]
SCszb = 17;			% Shouval [17]		Brad [17]

%}

SCdeltaT = 0.01;	% Shouval [.01]		Brad [.01]
SCbeta = 50;		% Shouval [60]		Brad [50]
SCtau = 1.8;		% Shouval [1.0]		Brad [1.8]
SCmu = 1/SCtau;		% Shouval [1/tau]	Brad [1/tau]
SCL = 1.2;			% Shouval [1.5]		Brad [1.2]
SCL2 = 1.1;			% Shouval [0.9]		Brad [1.1]
SCr = 15;			% Shouval [10]		Brad [15]
SCro = 0.90;		% Shouval [.95]		Brad [.90]
SCszi = 7;			% Shouval [8]		Brad [7]
SCszb = 17;			% Shouval [17]		Brad [17]

SC1szi = sap1;	SC2szi = sap2;
SC1beta = sap3;			SC2beta = sap6;		
SC1tau = sap4;			SC2tau = sap7;		
SC1mu = 1/SC1tau;			SC2mu = 1/SC2tau;	
SC1L = sap5;				SC2L = sap8;

SC1deltaT = 0.01;			SC2deltaT = 0.01;
SC1L2 = 1.1;				SC2L2 = 1.1;		
SC1r = 15;					SC2r = 15;			
SC1ro = 0.90;				SC2ro = 0.90;			
SC1szb = 17;				SC2szb = 17;

%===========================================%
%{
%==============%
%  S ORIGINAL
%==============%
time=Nsteps; % time is set in arbitrary units
delta_t=0.1; % time step, the smaller the more exact, and less parallel
num_epochs=floor(time/delta_t)+1; % 1001
initial_sz=7;
beta=60;
tau=1;
mu=1./tau; % (tau=1/mu), random degredation
r=10; % rate of transition of new receptor into unoccupied state
ro=0.95; %cellular concentration
ro2=ones(1,num_epochs)*0.00;
ro2(floor(15.5/delta_t):floor(16.5/delta_t))=0.02;
L1=1.5; % a good  number for repulsive lattice constant range ~ 1.2-2
L2=0.9;
sz_bound=17;
sz=sz_bound-2;
L1_array=L1*ones(1,num_epochs);
L2_array=L2*ones(1,num_epochs);
mu_array=mu*ones(1,num_epochs);
tau_2=1/4;
S=zeros(sz+2); % spin memory lattice
h_mask=[0 1 0; 1 0 1; 0 1 0];
center=floor(sz_bound/2)+1;
initial_down=center-floor(initial_sz/2);
initial_up=initial_down+initial_sz-1;
S(initial_down:initial_up,initial_down:initial_up)=1;
S1=S;
mem_num=zeros(1,num_epochs);
mem_num(1)=sum(S(:));
pop2_num=zeros(1,num_epochs);
% (S1pop2_num, h_mask, S1L1_array, S1L2_array, S1beta)
% (S1mu_array, S1delta_t, S1ro, S1r, S1ro2, S1, nn)
%}
%===========================================%
time=Nsteps;
h_mask=[0 1 0; 1 0 1; 0 1 0];
initial_sz=SCszi;
sz_bound=SCszb;
sz=sz_bound-2;
S=zeros(sz+2);
center=floor(sz_bound/2)+1;
initial_down=center-floor(initial_sz/2);
initial_up=initial_down+initial_sz-1;
S(initial_down:initial_up,initial_down:initial_up)=1;
%===========================================%
if runSAPPSD1 || runSAPPSD2 %   S1 PSD1
%===========================================%
sz1_bound=SC1szb;
sz1=sz1_bound-2;
S1=zeros(sz1+2);
center1=floor(sz1_bound/2)+1;
initial_down1=center1-floor(SC1szi/2);
initial_up1=initial_down1+SC1szi-1;
S1(initial_down1:initial_up1,initial_down1:initial_up1)=1;

S1delta_t=SC1deltaT; % smaller yields slower changes [0.1]
S1num_epochs=floor(time/S1delta_t)+1; % 1001
S1ro2=ones(1,S1num_epochs)*0.00;
S1ro2(floor(15.5/S1delta_t):floor(16.5/S1delta_t))=0.02;

S1beta=SC1beta;      % slope of on-rate 'h' probability insertion [60]
                % beta>100 low P() of psudopodia growth; beta<20 high P()

S1tau=SC1tau;      % off rate  random degredation [1] lower=fast degrade
S1mu=SC1mu;   % on rate   random degredation [1]
 
S1r=SC1r;          % receptor transition rate into unoccupied site [10]
S1ro=SC1ro;      % cellular concentration [.95]
 
% THIS MAKES THE CLUSTERS GROW OR SHRINK
S1L1=SC1L;       % repulsive lattice constant [1.2 - 2]
                % [3] MAKES IT SHRINK [1.2] MAKES IT GROW
S1L2=SC1L2;       % repulsive lattice constant [0.9] active during LTP window
 
S1L1_array=S1L1*ones(1,S1num_epochs);
S1L2_array=S1L2*ones(1,S1num_epochs);
S1mu_array=S1mu*ones(1,S1num_epochs);

S1s=S1;
S1mem_num=zeros(1,S1num_epochs);
S1mem_num(1)=sum(S1(:));
S1pop2_num=zeros(1,S1num_epochs);
%===========================================%
% end % if runSAPPSD1 %   S1 PSD1
%===========================================%

%===========================================%
% if runSAPPSD2 %   S2 PSD2
%===========================================%
sz2_bound=SC2szb;
sz2=sz2_bound-2;
S2=zeros(sz2+2);
center2=floor(sz2_bound/2)+1;
initial_down2=center2-floor(SC2szi/2);
initial_up2=initial_down2+SC2szi-1;
S2(initial_down2:initial_up2,initial_down2:initial_up2)=1;

S2delta_t=SC2deltaT; % smaller yields slower changes [0.1]
S2num_epochs=floor(time/S2delta_t)+1; % 1001
S2ro2=ones(1,S2num_epochs)*0.00;
S2ro2(floor(15.5/S2delta_t):floor(16.5/S2delta_t))=0.02;

S2beta=SC2beta;		% slope of on-rate 'h' probability insertion [60]
				% beta>100 low P() of psudopodia growth; beta<20 high P()

S2tau=SC2tau;		% off rate  random degredation [1] lower=fast degrade
S2mu=SC2mu;	% on rate	random degredation [1]

S2r=SC2r;			% receptor transition rate into unoccupied site [10]
S2ro=SC2ro;		% cellular concentration [.95]

% THIS MAKES THE CLUSTERS GROW OR SHRINK
S2L1=SC2L;		% repulsive lattice constant [1.2 - 2]
				% [3] MAKES IT SHRINK [1.2] MAKES IT GROW
S2L2=SC2L2;		% repulsive lattice constant [0.9]

S2L1_array=S2L1*ones(1,S2num_epochs);
S2L2_array=S2L2*ones(1,S2num_epochs);
S2mu_array=S2mu*ones(1,S2num_epochs);
S2s = S2;

S2mem_num=zeros(1,S2num_epochs);
S2mem_num(1)=sum(S2(:));
S2pop2_num=zeros(1,S2num_epochs);
%===========================================%
end % if runSAPPSD2 %   S2 PSD2
%===========================================%

if ~runSAPPSD1
S1=0;
end

if ~runSAPPSD2
S2=0;
end

%===========================================%
% CLUSTER SUMS FOR PSD POTENTIATION LEVEL
%-------------------------------------------%
S1sum = sum(S1(:)); S1sumOrig = sum(S1(:)); 
S2sum = sum(S2(:)); S2sumOrig = sum(S2(:));
%===========================================%
%-------------##########################------------------%



end





























