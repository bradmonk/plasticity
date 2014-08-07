function [varargout] = addaxis(varargin)


if nargin == 1
	
	if numel(varargin{1}) > 1
		axvals = varargin{1};
		XAxO = [axvals{1}(1) 0 0]; XAxT = [axvals{1}(2) 0 0];
		YAxO = [0 axvals{1}(3) 0]; YAxT = [0 axvals{1}(4) 0];
		ZAxO = [0 0 axvals{1}(5)]; ZAxT = [0 0 axvals{1}(6)];
	else
		axvals = varargin{1};
		XAxO = [-axvals 0 0]; XAxT = [axvals 0 0];
		YAxO = [0 -axvals 0]; YAxT = [0 axvals 0];
		ZAxO = [0 0 -axvals]; ZAxT = [0 0 axvals];
	end


else
	XAxO = [-15 0 0]; XAxT = [15 0 0];
	YAxO = [0 -15 0]; YAxT = [0 15 0];
	ZAxO = [0 0 -5]; ZAxT = [0 0 10];
end


[XAx YAx ZAx] = plot3prep({XAxT YAxT ZAxT},{XAxO YAxO ZAxO});
%-----------------
hAax = plot3(XAx,YAx,ZAx,'k','LineWidth',1,'MarkerSize',10);
set(hAax,{'LineWidth'},{1,1,1}'); hold on
varargout={hAax};
end