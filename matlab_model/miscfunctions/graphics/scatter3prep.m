function [varargout] = scatter3prep(varargin)
% Useage
% Ts = [x;y;z]
% [XMx YMx ZMx] = scatter3prep({Ts1 Ts2...})
%
% [XMx YMx ZMx] = scatter3prep({Ts1 Ts2...},{Mx})


%---------------
if nargin == 2
	Ts = varargin{1};
	Mx = varargin{3};
	Ts = cat(2,Ts,Mx);
	for nTs = 1:numel(Ts)
		TMx(nTs,:) = [Ts{nTs}(1);Ts{nTs}(2);Ts{nTs}(3)];
	end
%---------------
else
	Ts = varargin{1};
	for nTs = 1:numel(Ts)
		TMx(nTs,:) = [Ts{nTs}(1);Ts{nTs}(2);Ts{nTs}(3)];
	end
end
%---------------

XMx = TMx(:,1);
YMx = TMx(:,2);
ZMx = TMx(:,3);

varargout = {XMx,YMx,ZMx};
end