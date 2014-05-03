clc; close all; clear all;

% USER ENTERED VALUES
NSteps = 100;			% number of steps (loop time)
Ndots = 10;				% number of particles
DiffRate = 0.5;			% diffusion rate coefficient A

% BASE DIFFUSION RATES EQUATIONS
t = 1;						% time step (seconds)
dm = 2;                     % dimensions
D = DiffRate;				% Diffusion Rate A (D = L² / 2d*t)
k = sqrt(dm*D*t);			% stdev of D's step size distribution
MSD = 2*dm*D;				% theoretical mean squared displacement

% PREALLOCATE MATRICES FOR SPEED
XYL = zeros(2,Ndots);		% XY particle locations
XYS = zeros(2,Ndots);		% XY step sizes


% PERFORM BROWNIAN MOTION
for Nt = 1:NSteps 
	XYS = (k * randn(2,Ndots));	% generates step sizes
	XYL = XYL+XYS;				% adds step to location
	XL(Nt,:) = XYL(1,:);		% saves x location
	YL(Nt,:) = XYL(2,:);		% saves y location
end


% FORMAT MATRIX FOR JMP EXPORT
XLC = reshape(XL,[],1);
YLC = reshape(YL,[],1);
XYLC = [XLC YLC];
XYLDat = [XYLC ones(numel(XLC),1) ones(numel(YLC),1)];
TIME = repmat((1:100)',10,1);
XYLDat(:,3) = TIME;
PART = reshape(repmat((1:10),100,1),[],1);
XYLDat(:,4) = PART;
XYLD = cat(1,[0 0 0 0],XYLDat);
% xlswrite('DiffData',XYLD)
csvwrite('DiffData.csv',XYLD)




%{
Open(
	"/Users/bradleymonk/Documents/MatLab/BradsModel/miscfunctions/ASSORTED/DiffData.csv",
	columns(
		Column( "XL", Numeric, Continuous, Format( "Best", 10 ) ),
		Column( "YL", Numeric, Continuous, Format( "Best", 10 ) ),
		Column( "TIME", Numeric, Continuous, Format( "Best", 10 ) ),
		Column( "ID", Numeric, Continuous, Format( "Best", 10 ) )
	),
	Import Settings(
		End Of Line( CRLF, CR, LF ),
		End Of Field( Comma, CSV( 1 ) ),
		Strip Quotes( 0 ),
		Use Apostrophe as Quotation Mark( 0 ),
		Scan Whole File( 1 ),
		Treat empty columns as numeric( 0 ),
		CompressNumericColumns( 0 ),
		CompressCharacterColumns( 0 ),
		CompressAllowListCheck( 0 ),
		Labels( 0 ),
		Column Names Start( 0 ),
		Data Starts( 1 ),
		Lines To Read( "All" ),
		Year Rule( "20xx" )
	)
)
%}


