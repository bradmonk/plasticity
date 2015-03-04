function[]=mapdendrite(varargin)

if nargin > 0 
h=subplot(2,2,3);
set(h,'Position',[0.05 0.05 0.2 0.37])
fig1=imread('varargin{1}');
B=gray; mind=20;
B(1:mind,:)=B(1:mind,:)/B(mind,1);
colormap(1-B);
imagesc(fig1);
axis off;
else

h=subplot(2,2,3);
set(h,'Position',[0.05 0.05 0.2 0.37])
fig1=imread('DENDRITE.png');
B=gray; mind=20;
B(1:mind,:)=B(1:mind,:)/B(mind,1);
colormap(1-B);
imagesc(fig1);
axis off;
end;

