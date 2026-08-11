


OKPARAMS = LRPAR((AOK>0),3:6);
NOKPARAMS = LRPAR((AOK<1),3:6);

size(OKPARAMS) + size(NOKPARAMS)

[r,~,~] = find(OKPARAMS(:,1)==7);
QR7 = OKPARAMS(r,:);
[r,~,~] = find(OKPARAMS(:,1)==10.2);
QR10 = OKPARAMS(r,:);
[r,~,~] = find(OKPARAMS(:,1)==15);
QR15 = OKPARAMS(r,:);

QRMaxRows = max([size(QR7,1) size(QR10,1) size(QR10,1)]);
QRMinRows = min([size(QR7,1) size(QR10,1) size(QR10,1)]);

QRa = QR7;
QRa(end+1:QRMaxRows,:) = 0;
QRb = QR10;
QRb(end+1:QRMaxRows,:) = 0;
QRc = QR15;
QRc(end+1:QRMaxRows,:) = 0;


V=[];
V(:,:,1) = QRa(:,2:end);
V(:,:,2) = QRb(:,2:end);
V(:,:,3) = QRc(:,2:end);

Xs = [QRa(1) QRb(1) QRc(1)];
Sx = repmat(Xs,[QRMaxRows 1]);
Ys = linspace(7,15,QRMaxRows)';
Sy = repmat(Ys,[1 3]);
Zs = ones(size(Sx));
Sz = Zs + min(min(OKPARAMS));

xx = 1:numel(V(:,1,1));
yy = 1:numel(V(1,:,1));
zz = 1:numel(V(1,1,:));

slice(V,xx,yy,zz)
surf(V)




Xs = OKPARAMS(:,2); 
Ys = OKPARAMS(:,3);
Zs = OKPARAMS(:,4);

Xm = meshgrid(Xs);
Ym = meshgrid(Ys);

surf(Xs,Ys,Zs)

















