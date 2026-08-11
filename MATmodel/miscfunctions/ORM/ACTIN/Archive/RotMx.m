function [] = RotMx()
clc; close all; clear all;

% M=M3d(AngleDegrees, FirstPoint-SecondPoint, SecondPoint); %transformation matrix
M=M3d([45 45 70]', [5 5 7]', [10 10 3]');

% NewPoint=M(1:3,:) * [ThirdPoint(:);1];

end



function M=M3d(deg,u,x0)
%Generate roto-translation matrix for the rotation around an arbitrary line in 3D.
%
% M=M3d(deg,u,x0)
%
%in:
%
% deg: The counter-clockwise rotation about the line in degrees.
% u,x0: 3D vectors specifying the line in parametric form x(t)=x0+t*u 
% Default for x0=0 (pure rotation).
%out:
%
% M: A 4x4 homogenous coordinate transform matrix representing
% the roto-translation. 
%
%See also: Rx,Ry,Rz, R2d, R3d, M2d


if nargin<2, x0=[0;0;0]; end

x0=x0(:); u=u(:)/norm(u);

AxisShift=x0-(x0.'*u).*u;




Mshift=mkaff(eye(3),-AxisShift);

Mroto=mkaff(R3d(deg,u));

M=inv(Mshift)*Mroto*Mshift;

end


function M=mkaff(R,t)

if nargin<2, t=[0;0;0]; end

nn=size(R,1);


M=eye(nn+1);

M(1:end-1,1:end-1)=R;
M(1:end-1,end)=t(:);

end


function R=R3d(deg,u)
%R3D - 3D Rotation matrix counter-clockwise about an axis.
%
%R=R3d(deg,axis)
%
% deg: The counter-clockwise rotation about the axis in degrees.
%
%See also Rx,Ry,Rz,R3d,M2d,M3d

R=eye(3);
u=u(:)/norm(u);
x=deg; %abbreviation

for ii=1:3
   
    v=R(:,ii);
    
    R(:,ii)=v*cosd(x) + cross(u,v)*sind(x) + (u.'*v)*(1-cosd(x))*u;
      %Rodrigues' formula
      
end

end