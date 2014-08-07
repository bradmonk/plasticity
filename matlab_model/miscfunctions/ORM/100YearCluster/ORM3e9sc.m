

profile on;
ORM3e9()
profile viewer;






%{
hMon = hkMask*Bon;
LBon = Lon*Bon;

hMoff = hkMask*Boff;
LBoff = Loff*Boff;


Pmx = rand(size(S));
Soc = (S>0);

hon = -convn(Soc,hMon,'same')+LBon;
hoff = convn(Soc,hMoff,'same')-LBoff;

Pkon = ~Soc .* (dTRon ./ (1+exp(hon)) );
Pkoff = Soc .* (dTRoff ./ (1+exp(hoff)) );

Son = (Pkon>Pmx);
Soff = (Pkoff>Pmx);

S = (Soc-Soff) + Son;
%}

%{


hk = convn(Soc,hkMask,'same');
Lhn = (hk-Lon) .* (-Bon);
Pkn = ~Soc .* ( dTRon ./ (1+exp(Lhn)) );
Lhf = ((-hk)+Loff) .* (-Boff);
Pkf = Soc .* ( dTRoff ./ (1+exp(Lhf)) );
Sn = (Pkn>Pmx);
Sf = (Pkf>Pmx);

%--
Pmx = rand(size(S));
Soc = (S>0);
hk = convn(Soc,hkMask,'same');


Lhon = (hk-Lon) .* (-Bon);
Pkon = ~Soc .* ( dTRon ./ (1+exp(Lhon)) );

Lhoff = ((-hk)+Loff) .* (-Boff);
Pkoff = Soc .* ( dTRoff ./ (1+exp(Lhoff)) );

Son = (Pkon>Pmx);
Soff = (Pkoff>Pmx);

S = (Soc-Soff) + Son;
%}

%{
Pmx = rand(size(S));
Soc = (S>0);
Sno = ~Soc;
hk = convn(Soc,hkMask,'same');

Lhon = (hk-Lon) .* (-Bon);
Pon = 1 ./ (1+exp(Lhon));
Pkon = Sno .* ( dT * Ron * Pon );
Son = (Pkon>Pmx);

Lhoff = ((-hk)+Loff) .* (-Boff);
Poff = 1 ./ (1+exp(Lhoff));
Pkoff = Soc .* ( dT * Roff * Poff );
Soff = (Pkoff>Pmx);

S = (Soc-Soff) + Son;
%}



%{

%}



%{

%}


%{

%}


%{

%}


%{

%}

























