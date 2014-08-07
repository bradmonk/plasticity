function [] = anotext(varargin)

tbgc = [.91 .91 .91];
nt = nargin;

% inputname(1)

for an = 1:nt
val=[];pos=[];txt=[];

if numel(varargin{an}) >= 1
val = varargin{an}{1};
if numel(val) == 3
val=[val(1) val(2) val(3)];
else
val=[0 0 0]+an;
end
else
val=[0 0 0]+an;
end

if numel(varargin{an}) >= 2
txt = varargin{an}{2};
txt = strcat('\leftarrow ',txt);
else
txt='\bullet';
end

if numel(varargin{an}) >= 3
pos = varargin{an}{3};
	if numel(pos) == 3
	pos=[pos(1) pos(2) pos(3)];
	else
	pos=val;
	end
else
pos=val;
end

% hT = text(pos(1),pos(2),pos(3),...
% strcat('\bullet ', txt,'(', sprintf('%.2f|',val'),')'),...
% 'FontSize',12,'HorizontalAlignment','left');

hT = text(pos(1),pos(2),pos(3),...
strcat('\bullet ', txt,'(', sprintf('%.2f|',val'),')'),...
'FontSize',12,'HorizontalAlignment','left',...
'BackgroundColor',tbgc);

end

% hBGC = findobj('BackgroundColor',tbgc);
alpha(.05);

end