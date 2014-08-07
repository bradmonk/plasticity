
dT = .01;
Lon = 2;
Loff = 2;
Bon = 12;
Boff = 1;
Ron = 16;
Roff = 4;


S=padarray(ones(5),[2 2], 0);
hkMask=[0 1 0; 1 0 1; 0 1 0];

dT = .01;
Lon = 3;
Loff = 2;
Bon = 2;
Boff = 1;
Ron = 1;
Roff = 1;

Pmx = rand(size(S));
Soc = (S>0);
Sno = ~Soc;
hk = convn(Soc,hkMask,'same');

Pon = 1 ./ (1+exp((hk-Lon).*(-Bon)));
Pkon = Sno .* ( Ron * dT * Pon );
Son = (Pkon>Pmx);

Poff = 1 ./ (1+exp(((-hk)+Loff).*(-Boff)));
Pkoff = Soc .* ( Roff * dT * Poff );
Soff = (Pkoff>Pmx);


hLB = exp((hk-Lon).*(-Bon));

xx = [0 1 2 3 4];
Sz = [S(3,1) S(3,2) S(3,3) S(4,3) S(4,4)];
hkz = [hk(3,1) hk(3,2) hk(3,3) hk(4,3) hk(4,4)];
hLBz = [hLB(3,1) hLB(3,2) hLB(3,3) hLB(4,3) hLB(4,4)];


cc = 0;
for nn = 1:3
cc = cc+.1;
S=padarray(ones(5),[2 2], 0);
hkMask=[0 1 0; 1 0 1; 0 1 0];

dT = .01;
Lon = 3;
Bon = nn/2;
Ron = 1;

Soc = (S>0); Sno = ~Soc; hk = convn(Soc,hkMask,'same');

hLB = exp((hk-Lon).*(-Bon));
Pon = 1 ./ (1+exp((hk-Lon).*(-Bon)));

xx = [0 1 2 3 4];
hLBz = [hLB(3,1) hLB(3,2) hLB(3,3) hLB(4,3) hLB(4,4)];
Po = [Pon(3,1) Pon(3,2) Pon(3,3) Pon(4,3) Pon(4,4)];
LB{nn} = hLBz;

figure(1)
% plot(xx,Po,':','Color',[(.9-cc*1.5) (.3) (.5+cc)],'LineWidth',2)
% plot(xx,Po,'Color',[(.1+cc) (.1) (.9)],'LineWidth',2)
% hold on

% subplot(2,1,1),
% plot(xx,hLBz,'Color',[(.9-cc*1.5) (.1+cc) (.5+cc)],'LineWidth',2)
% hold on
% subplot(2,1,2),
% plot(xx,Po,':','Color',[(.9-cc*1.5) (.3) (.5+cc)],'LineWidth',2)
% plot(xx,Po,'Color',[(.1+cc) (.1) (.9)],'LineWidth',2)
% hold on

pause(.5)
end


















cc = 0;
for nn = 1:3
S=padarray(ones(5),[2 2], 0);
hkMask=[0 1 0; 1 0 1; 0 1 0];

dT = .01;
Lon = 3;
Bon = nn/2;
Ron = 1;

Lon1=1;
Lon3=3;

Soc = (S>0); Sno = ~Soc; hk = convn(Soc,hkMask,'same');

Pon1 = 1 ./ (1+exp((hk-Lon1).*(-Bon)));
Pon3 = 1 ./ (1+exp((hk-Lon3).*(-Bon)));

Po1 = [Pon1(3,1) Pon1(3,2) Pon1(3,3) Pon1(4,3) Pon1(4,4)];
Po3 = [Pon3(3,1) Pon3(3,2) Pon3(3,3) Pon3(4,3) Pon3(4,4)];

xx = [0 1 2 3 4];
figure(1)
plot(xx,Po1,'-.','Color',[(1-cc) (0) (1)],'LineWidth',5)
hold on
plot(xx,Po3,'Color',[(1-cc) (0) (1)],'LineWidth',5)
hold on

pause(.5)
cc = cc+.3;
end








clc; close all; clear all;
cc = 0;
for nn = 1:3
S=padarray(ones(5),[2 2], 0);
hkMask=[0 1 0; 1 0 1; 0 1 0];

dT = .1;
Lon = 2;
Ron = nn*2;

Bon1 = 1;
Bon2 = 3;


Soc = (S>0); Sno = ~Soc; hk = convn(Soc,hkMask,'same');

Pon1 = 1 ./ (1+exp((hk-Lon).*(-Bon1)));
Pon2 = 1 ./ (1+exp((hk-Lon).*(-Bon2)));

Pkon1 = ( Ron * dT * Pon1 );
Pkon2 = ( Ron * dT * Pon2 );

Po1 = [Pkon1(3,1) Pkon1(3,2) Pkon1(3,3) Pkon1(4,3) Pkon1(4,4)];
Po2 = [Pkon2(3,1) Pkon2(3,2) Pkon2(3,3) Pkon2(4,3) Pkon2(4,4)];

xx = [0 1 2 3 4];
figure(1)
hPo1 = plot(xx,Po1,'-.','Color',[(1-cc) (0) (1)],'LineWidth',5);
hold on
plot(xx,Po2,'Color',[(1-cc) (0) (1)],'LineWidth',5)
hold on

pause(.5)
cc = cc+.3;
end
% legend(hPo1,'Bon1 Ron1','Bon1 Ron2','Bon1 Ron3')



% [1 0 0] red
% 
% [1 0 1] magenta 
% 
% [0 0 1] blue
% 
% [0 1 0] green

