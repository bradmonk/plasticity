


time=800;		% (200) time is set in arbitrary units
delta_t=0.1;	% (0.025) time step, the smaller the more exact, and less parallel
num_epochs=floor(time/delta_t)+1; % 8001


beta=50;		% (60) insertion probability slope (large=stable, small=growth)
tau=1;			% (1) random degredation (big=5=filled cluster)
mu=1./tau;		% (mu=1/tau) internalization rate
r=10;			% (10) rate of transition of new receptor into unoccupied state
ro=0.85;		% (.95) cellular concentration
L1=1.5;			% (1.5) lattice repulsion constant

h_mask=[0 1 0; 1 0 1; 0 1 0];

Csz = 7;
Bsz = 3;
S=ones(Csz);
S=padarray(S,[Bsz Bsz]);
S1=S;


Tsz = numel(S);			% Total matrix space
Ssz = numel(find(S>0));	% Filled matrix space
Nsz = Tsz-Ssz;			% Empty matrix space
Psz = round(Ssz/Tsz*100);	% Percent filled space
L1x = linspace(0.8,3.5);



mem_num=zeros(1,num_epochs);
mem_num(1)=sum(S(:));
pop2_num=zeros(1,num_epochs);
figure;

for nn=2:num_epochs
	S_occ=(S>0);  
	mem_num(nn)=sum(S_occ(:));
	h=convn(S_occ,h_mask,'same');
	
	Tsz = numel(S);			% Total matrix space
	Ssz = numel(find(S>0));	% Filled matrix space
	Psz = round(Ssz/Tsz*100);	% Percent filled space
	L1 = L1x(Psz);

	dE1=L1-h;					% energy diff
	P1=1./(1+exp(dE1*beta));	% conditional exocytosis

	Pen=S_occ.*(1-exp(-mu*delta_t)); % large delta_t approx
	P_rand=rand(size(S));
	S_en=(Pen>P_rand);

	P_ex1=(1-S_occ).*(ro*r*delta_t*P1);
	S_ex=(P_ex1>P_rand);
	S=(S_occ-S_en).*S+S_ex;
 
 
 if mod(nn,10)==0
 imagesc(S)
 colormap bone;
 drawnow
 end
 
end % for nn


time_line=[0:delta_t:delta_t*(num_epochs-1)];
figure
plot(time_line,mem_num);
hold on
plot(time_line,pop2_num,'r');
  
  