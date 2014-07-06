
kk=0;a=0.1;N=300;zeta=10;DD=0.0001;rrr=0.01;ar=0.5; 
% parameters in node strength; aver filam length; init number of asters; rate of nucl
h=0.25; NN=300/h; % time steps; keep rate*N < 1
x=sqrt(N)*rand(N,1);y=sqrt(N)*rand(N,1); % coordinates of moving asters
dispx=zeros(size(x));dispy=zeros(size(x));
x1=zeros(size(1:floor(sqrt(N))));y1=sqrt(N)*rand(1,length(x1)); % coord of fixed asters at the edge
x2=sqrt(N)*ones(size(1:floor(sqrt(N))));y2=sqrt(N)*rand(1,length(x2)); % coord of fixed asters at the edge
y3=zeros(size(1:floor(sqrt(N))));x3=sqrt(N)*rand(1,length(y3)); % coord of fixed asters at the edge
y4=sqrt(N)*ones(size(1:floor(sqrt(N))));x4=sqrt(N)*rand(1,length(y4)); % coord of fixed asters at the edge
xx=[x1, x2, x3, x4];yy=[y1, y2, y3, y4]; % coord of all fixed asters
zs=ones(size(xx)); z=ones(size(x)); % strength of aster
dob=zeros(size(x)); % date of birth of asters
Vx=zeros(size(x));Vy=zeros(size(x)); % vel of asters
Int=[];SSS=[];XXZ=[];VV=[];Oold=[];Mmsd=[]; % numbers of interact, mergers, nodes, aver vel
 
for g=1:NN
    Epp=0;SS=0;V=0;
j=1;while j<length(x)+1
    for s=1:length(xx) % summing interactions for a given node with fixed nodes
        rr=sqrt((x(j)-xx(s)).^2+(y(j)-yy(s)).^2)+0.01; % distance between nodes
        ep=0.5*(sign(exp(-rr/a)-rand)+1)*sign(rand-ar); % length of actin filaments is distributed exponentially;
        % plus a random factor; ep is variable showing if two nodes interact at all
        Vx(j)=Vx(j)-ep*(x(j)-xx(s))/rr; % adding x-component of interaction
        Vy(j)=Vy(j)-ep*(y(j)-yy(s))/rr; % adding y-component of interaction
    end
    k=j;while k<length(x)+1 % summing interactions for a given node with other moving nodes
        r=sqrt((x(j)-x(k)).^2+(y(j)-y(k)).^2)+0.01; % distance between nodes
        %ep=0.5*(sign(-((r.^2)/(r.^2+a^2))+rand)+1); tr=z(j)*z(k)/(z(j)*z(k)+kk); Epp=Epp+ep;
        ep=0.5*(sign(exp(-r/a)-rand)+1)*sign(rand-ar); tr=z(j)*z(k)/(z(j)*z(k)+kk); Epp=Epp+ep; %length of filaments,node strength
        Vx(j)=Vx(j)-tr*ep*(x(j)-x(k))/r;Vx(k)=Vx(k)+tr*ep*(x(j)-x(k))/r;% adding x-comp of interaction
        Vy(j)=Vy(j)-tr*ep*(y(j)-y(k))/r;Vy(k)=Vy(k)+tr*ep*(y(j)-y(k))/r;% adding y-comp of interaction
        k=k+1; z(j)=z(j)-h*0.05*(z(j)-1); % evolution of node strength
    end
    dddx=h*Vx(j)/(zeta*z(j))+sqrt(2*h*DD)*randn;x(j)=x(j)+dddx;dispx(j)=dispx(j)+dddx;
    dddy=h*Vy(j)/(zeta*z(j))+sqrt(2*h*DD)*randn;y(j)=y(j)+dddy;dispy(j)=dispy(j)+dddy;
    V=V+sqrt(Vx(j)^2+Vy(j)^2)/(zeta*z(j)); dob(j)=dob(j)+h; j=j+1;
end
    j=1;while j<length(x)+1 
        c=x(j);b=y(j);zzz=sum(z);ddx=dispx(j);ddy=dispy(j);% updating node position 
        s=find(abs(x-c)+abs(y-b)>rrr);z=z(s);zz=sum(z); % finding all nodes that are not too close to the one
        SS=SS+length(x)-length(s)-1; if length(x)>length(s)+1 dob=[dob(s); 0]; else dob=dob; end
        x=x(s);y=y(s);dispx=dispx(s);dispy=dispy(s);
        x=[x; c];y=[y; b]; dispx=[dispx; ddx];dispy=[dispy; ddy];z=[z; (zzz-zz)]; 
        % merging nodes that are too close, pooling their strength; updating age
    j=j+1; end
    if zzz<N x=[x; sqrt(N)*rand];y=[y; sqrt(N)*rand];z=[z; 1];dispx=[dispx; 0];dispy=[dispy; 0];dob=[dob; 0]; else 1==1; end 
% nucleation of new nodes
    Vx=zeros(size(x));Vy=zeros(size(x)); % resetting velocities for the next loop
    XZ=length(x);V=V/length(x);Int=[Int Epp];SSS=[SSS SS];XXZ=[XXZ XZ];VV=[VV V];
    sss=find(dob==h*g);Old=length(sss); Oold=[Oold Old];dix=dispx(sss);diy=dispy(sss);
    mssd=sum((dispx.^2+dispy.^2))/length(sss);Mmsd=[Mmsd mssd];
    % numbers of interact, mergers, nodes, aver vel
end
 
% plot(h*(1:NN),XXZ/N,'k',h*(1:NN),SSS,'m',h*(1:NN),Int./XXZ,'r');   plot(h*(1:NN),Oold,h*(1:NN),Mmsd)
% plot(log(1:NN),log(Mmsd)-log(Mmsd(1)),log(1:NN),0.5*log(2*DD*h*(1:NN))-0.5*log(2*DD*h))
