function [] = ExampleMainFunction()
clc; close all; clear all;


%---------------------------------------------%
% a function of of seven variables that yields the result of
% rotating the point (x,y,z) about the axis ?u,v,w? by the angle ?.
Fxut = @(x,y,z,u,v,w,tta) ...
	([(u*(u*x+v*y+w*z)*(1-cos(tta))+(u^2+v^2+w^2)*...
	x*cos(tta)+sqrt(u^2+v^2+w^2)*(-w*y+v*z)*sin(tta))/... 
	(u^2 + v^2 + w^2);
	(v*(u*x+v*y+w*z)*(1-cos(tta))+(u^2+v^2+w^2)*...
	y*cos(tta)+sqrt(u^2+v^2+w^2)*(w*x-u*z)*sin(tta))/... 
	(u^2+v^2+w^2);
	(w*(u*x+v*y+w*z)*(1-cos(tta))+(u^2+v^2+w^2)*...
	z*cos(tta)+sqrt(u^2+v^2+w^2)*(-v*x+u*y)*sin(tta))/... 
	(u^2+v^2+w^2)]);
% MFxut = Fxut(1,2,3,4,5,6,.8)
%---------------------------------------------%



%---------------------------------------------%
% We will define an arbitrary line by a point the line goes through and a 
% direction vector. If the axis of rotation is given by two points P1 = (a,b,c) 
% and P2 = (d,e,f), then a direction vector can be obtained by 
% ?u,v,w? = ?d-a,e-b,f-c?
% We can now write a transformation for the rotation of a point about this line.
% 
% The matrix for rotation about an arbitrary line
% If we multiply this times ?x,y,z? we can obtain a function of of ten variables 
% that yields the result of rotating the point (x,y,z) about the 
% line through (a,b,c) with direction vector ?u,v,w? by the angle ?.
Fxaut = @(x,y,z,a,b,c,u,v,w,tta) ... 
	([(((a*(v^2+w^2)-u*(b*v+c*w-u*x-v*y-w*z))*(1-cos(tta))+(u^2+v^2+w^2)*x*cos(tta)+...
	sqrt((u^2+v^2+w^2))*(-c*v+v*w-w*y+v*z)*sin(tta))/(u^2+v^2+w^2));
	(((b*(u^2+w^2)-v*(a*u+c*w-u*x-v*y-w*z))*(1-cos(tta))+(u^2+v^2+w^2)*y*cos(tta)+...
	sqrt((u^2+v^2+w^2))*(c*u-a*w+w*x-u*z)*sin(tta))/(u^2+v^2+w^2));
	(((c*(u^2+v^2)-w*(a*u+c*w-u*x-v*y-w*z))*(1-cos(tta))+(u^2+v^2+w^2)*z*cos(tta)+...
	sqrt((u^2+v^2+w^2))*(-b*u+a*v-v*x+u*y)*sin(tta))/(u^2+v^2+w^2))]);
% MFxaut = Fxaut(5,0,0,1,1,1,1,1,1,.8)
%---------------------------------------------%




xyz = [3;3;3];
abc = [5;0;0];
def = [1;1;1];
uvw = abc-def;
% uvw = uvw ./ sqrt(sum(uvw.^2));

abcT = abc;
abcO = def;

x=xyz(1);y=xyz(2);z=xyz(3);
a=abc(1);b=abc(2);c=abc(3);
u=uvw(1);v=uvw(2);w=uvw(3);
UVW0 = repmat(uvw,1,50);

for ang = 1:50
	tta = rand*2*pi;
	RMX1(:,ang) = Fxut(x,y,z,u,v,w,tta);
end

Xop1 = [UVW0(1,:); RMX1(1,:)];
Yop1 = [UVW0(2,:); RMX1(2,:)];
Zop1 = [UVW0(3,:); RMX1(3,:)];

for ang = 1:50
	tta = rand*2*pi;
	RMX2(:,ang) = Fxaut(x,y,z,a,b,c,u,v,w,tta);
end

Xop2 = [UVW0(1,:); RMX2(1,:)];
Yop2 = [UVW0(2,:); RMX2(2,:)];
Zop2 = [UVW0(3,:); RMX2(3,:)];


% figure(1)
% plot3(Xop1,Yop1,Zop1)
% hold on
% plot3(Xop2,Yop2,Zop2)

figure(2)
% scatter3(Xop1(2,:),Yop1(2,:),Zop1(2,:),'.r')
% hold on
scatter3(Xop2(2,:),Yop2(2,:),Zop2(2,:),'.b')
hold on
scatter3(xyz(1),xyz(2),xyz(3),'o')
hold on
scatter3(abc(1),abc(2),abc(3),'d')
hold on
scatter3(uvw(1),uvw(2),uvw(3),'x')
axis([-15 15 -15 15 -15 15])
hold on

sFilO = [-15 0 0; 0 -15 0; 0 0 -5];
sFilT = [15 0 0; 0 15 0; 0 0 15];
sFilO = cat(2,abcO,sFilO);
sFilT = cat(2,abcT,sFilT);
Xop0 = [sFilO(1,:); sFilT(1,:)];
Yop0 = [sFilO(2,:); sFilT(2,:)];
Zop0 = [sFilO(3,:); sFilT(3,:)];
plot3(Xop0,Yop0,Zop0)
axis vis3d
%---------------------------------------------%

keyboard

end






