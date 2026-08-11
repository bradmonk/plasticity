function [XMx YMx ZMx] = plot3prep(varargin)
% Useage
% [XMx YMx ZMx] = plot3prep({Ts},{Os})
% [XMx YMx ZMx] = plot3prep({T1 T2 T3 T4},{O1 O2 O3 O4});
%
% [XMx YMx ZMx] = plot3prep({Ts},{Os},{Mx})
% [XMx YMx ZMx] = plot3prep({T1 T2},{O1 O2},{Mx});


%---------------
if nargin == 3
	Ts = varargin{1};
	Os = varargin{2};
	Mx = varargin{3};
	Ts = cat(2,Ts,Mx);
	Os = cat(2,Os,Mx);
	for nTs = 1:numel(Ts)
		XMx(:,nTs) = [Os{nTs}(1); Ts{nTs}(1)];
		YMx(:,nTs) = [Os{nTs}(2); Ts{nTs}(2)];
		ZMx(:,nTs) = [Os{nTs}(3); Ts{nTs}(3)];
	end
%---------------
else
	Ts = varargin{1};
	Os = varargin{2};
	for nTs = 1:numel(Ts)
		XMx(:,nTs) = [Os{nTs}(1); Ts{nTs}(1)];
		YMx(:,nTs) = [Os{nTs}(2); Ts{nTs}(2)];
		ZMx(:,nTs) = [Os{nTs}(3); Ts{nTs}(3)];
	end
end
%---------------


% varargout = {XMx,YMx,ZMx};
end