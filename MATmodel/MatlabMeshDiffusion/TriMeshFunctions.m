function varargout = TriMeshFunctions(varargin)


[argu, itm] = deal(varargin{:});

if strcmpi(argu, 'get_spine')
    varargout = {get_spine(itm)}; return
end


if strcmpi(argu, 'get_dendrite')
    varargout = {get_dendrite(itm)}; return
end


%% BUILD SPINE SKELETON

function SPINE_ALL = get_spine(itm)

  SCALE_FACTOR = 100;    % Mesh has issues below a certain scale

  if (itm == 1)
    % %---------------- SPINE OBJECT SHAPE A -----------------%
    TOTAL_LENGTH = 1.4;     TL=TOTAL_LENGTH; % total length of spine

    lt = @(x,a,b,c,d) (SCALE_FACTOR*(TL-(c*(1-(x-a)/(b-a)) + d*((x-a)/(b-a)))));

    % Z-AXIS VALUES FOR HIGHEST POINT OF SPINE OBJECTS
    Z_SYNTOP = lt(   0.00   ,0,1,0,TL);  % SYNAPSE_SURFACE
    Z_PERPSD = lt(   0.08   ,0,1,0,TL);  % PERISYNAPTIC_AREA
    Z_HEDTOP = lt(   0.20   ,0,1,0,TL);  % SPINE_HEAD_UPPER
    Z_HEDMID = lt(   0.30   ,0,1,0,TL);  % SPINE_HEAD_MIDDLE
    Z_HEDLOW = lt(   0.40   ,0,1,0,TL);  % SPINE_HEAD_LOWER
    Z_SPNECK = lt(   0.52   ,0,1,0,TL);  % SPINE_NECK
    Z_NKLINK = lt(   0.91   ,0,1,0,TL);  % NECK_LINKER
    Z_SHLINK = lt(   1.00   ,0,1,0,TL);  % SHAFT_LINKER
    
    % RADIUS AT TOP & BOTTOM OF SPINE OBJECTS
    R_SYNTOP = SCALE_FACTOR * 0.22 ;    % SURF SYNAPSE_SURFACE
    R_tPERPSD = SCALE_FACTOR * 0.30 ;   % TOP PERISYNAPTIC_AREA
    R_bPERSYN = SCALE_FACTOR * 0.35;   % BOT PERISYNAPTIC_AREA
    R_tHEDTOP = R_bPERSYN          ;  % TOP SPINE_HEAD_UPPER
    R_bHEDTOP = SCALE_FACTOR * 0.30;   % BOT SPINE_HEAD_UPPER
    R_tHEDMID = R_bHEDTOP          ;   % TOP SPINE_HEAD_MIDDLE
    R_bHEDMID = SCALE_FACTOR * 0.18;   % BOT SPINE_HEAD_MIDDLE
    R_tHEDLOW = R_bHEDMID          ;   % TOP SPINE_HEAD_LOWER
    R_bHEDLOW = SCALE_FACTOR * 0.12;   % BOT SPINE_HEAD_LOWER
    R_tSPNECK = R_bHEDLOW          ;   % TOP SPINE_NECK
    R_bSPNECK = SCALE_FACTOR * 0.10;   % BOT SPINE_NECK
    R_tNKLINK = R_bSPNECK          ;   % TOP NECK_LINKER
    R_bNKLINK = SCALE_FACTOR * 0.16;   % BOT NECK_LINKER
    R_tSHLINK = R_bNKLINK          ;   % TOP SHAFT_LINKER
    R_bSHLINK = R_tSHLINK          ;   % BOT SHAFT_LINKER
    % (LAST TWO MUST BE THE SAME VALUE)

            [xp, yp] = circus(R_SYNTOP, 1 , 1);
            zp = zeros(size(xp)) + Z_SYNTOP;
        SPINE_PSD = [xp; yp; zp]';

            [xp, yp] = circus(R_tPERPSD, 1 , 1);
            zp = zeros(size(xp)) + Z_PERPSD;
        SPINE_PERPSD = [xp; yp; zp]';

            [xp, yp] = circus(R_bPERSYN, 1 , 1);
            zp = zeros(size(xp)) + Z_HEDTOP;
        SPINE_SYNAPSE = [xp; yp; zp]';

            [xp, yp] = circus(R_tHEDMID, 1 , 1);
            zp = zeros(size(xp)) + Z_HEDMID;
        SPINE_HEAD_UPPER = [xp; yp; zp]';

            [xp, yp] = circus(R_tHEDLOW, 1 , 1);
            zp = zeros(size(xp)) + Z_HEDLOW;
        SPINE_HEAD_LOWER = [xp; yp; zp]';

            [xp, yp] = circus(R_tSPNECK, 1 , 1);
            zp = zeros(size(xp)) + Z_SPNECK;
        SPINE_NECK = [xp; yp; zp]';

            [xp, yp] = circus(R_tNKLINK, 1 , 1);
            zp = zeros(size(xp)) + Z_NKLINK;
        NECK_LINKER = [xp; yp; zp]';

            [xp, yp] = circus(R_tSHLINK, 1 , 1);
            zp = zeros(size(xp)) + Z_SHLINK;
        SHAFT_LINKER = [xp; yp; zp]';


    SPINE_ALL = [SPINE_PSD; SPINE_PERPSD; SPINE_SYNAPSE; SPINE_HEAD_UPPER; ... 
                 SPINE_HEAD_LOWER; SPINE_NECK; NECK_LINKER; SHAFT_LINKER];
  else

    % %---------------- SPINE OBJECT SHAPE B -----------------%
    TOTAL_LENGTH = 1.0;     TL=TOTAL_LENGTH; % total length of spine

    lt = @(x,a,b,c,d) (SCALE_FACTOR*(TL-(c*(1-(x-a)/(b-a)) + d*((x-a)/(b-a)))));

    % Z-AXIS VALUES FOR HIGHEST POINT OF SPINE OBJECTS
    Z_SYNTOP = lt(   0.00   ,0,1,0,TL);  % SYNAPSE_SURFACE
    Z_PERPSD = lt(   0.10   ,0,1,0,TL);  % PERISYNAPTIC_AREA
    Z_HEDTOP = lt(   0.25   ,0,1,0,TL);  % SPINE_HEAD_UPPER
    Z_HEDMID = lt(   0.35   ,0,1,0,TL);  % SPINE_HEAD_MIDDLE
    Z_HEDLOW = lt(   0.45   ,0,1,0,TL);  % SPINE_HEAD_LOWER
    Z_SPNECK = lt(   0.60   ,0,1,0,TL);  % SPINE_NECK
    Z_NKLINK = lt(   0.91   ,0,1,0,TL);  % NECK_LINKER
    Z_SHLINK = lt(   1.00   ,0,1,0,TL);  % SHAFT_LINKER

    % RADIUS AT TOP & BOTTOM OF SPINE OBJECTS
    R_SYNTOP = SCALE_FACTOR * 0.22 ;    % SURF SYNAPSE_SURFACE
    R_tPERPSD = SCALE_FACTOR * 0.32;   % TOP PERISYNAPTIC_AREA
    R_bPERSYN = SCALE_FACTOR * 0.38;   % BOT PERISYNAPTIC_AREA
    R_tHEDTOP = R_bPERSYN          ;  % TOP SPINE_HEAD_UPPER
    R_bHEDTOP = SCALE_FACTOR * 0.35;   % BOT SPINE_HEAD_UPPER
    R_tHEDMID = R_bHEDTOP          ;   % TOP SPINE_HEAD_MIDDLE
    R_bHEDMID = SCALE_FACTOR * 0.28;   % BOT SPINE_HEAD_MIDDLE
    R_tHEDLOW = R_bHEDMID          ;   % TOP SPINE_HEAD_LOWER
    R_bHEDLOW = SCALE_FACTOR * 0.17;   % BOT SPINE_HEAD_LOWER
    R_tSPNECK = R_bHEDLOW          ;   % TOP SPINE_NECK
    R_bSPNECK = SCALE_FACTOR * 0.11;   % BOT SPINE_NECK
    R_tNKLINK = R_bSPNECK          ;   % TOP NECK_LINKER
    R_bNKLINK = SCALE_FACTOR * 0.14;   % BOT NECK_LINKER
    R_tSHLINK = R_bNKLINK          ;   % TOP SHAFT_LINKER
    R_bSHLINK = R_tSHLINK          ;   % BOT SHAFT_LINKER
    % (LAST TWO MUST BE THE SAME VALUE)

            [xp, yp] = circus(R_SYNTOP, 1 , 1);
            zp = zeros(size(xp)) + Z_SYNTOP;
        SPINE_PSD = [xp; yp; zp]';

            [xp, yp] = circus(R_tPERPSD, 1 , 1);
            zp = zeros(size(xp)) + Z_PERPSD;
        SPINE_PERPSD = [xp; yp; zp]';

            [xp, yp] = circus(R_bPERSYN, 1 , 1);
            zp = zeros(size(xp)) + Z_HEDTOP;
        SPINE_SYNAPSE = [xp; yp; zp]';

            [xp, yp] = circus(R_tHEDMID, 1 , 1);
            zp = zeros(size(xp)) + Z_HEDMID;
        SPINE_HEAD_UPPER = [xp; yp; zp]';

            [xp, yp] = circus(R_tHEDLOW, 1 , 1);
            zp = zeros(size(xp)) + Z_HEDLOW;
        SPINE_HEAD_LOWER = [xp; yp; zp]';

            [xp, yp] = circus(R_tSPNECK, 1 , 1);
            zp = zeros(size(xp)) + Z_SPNECK;
        SPINE_NECK = [xp; yp; zp]';

            [xp, yp] = circus(R_tNKLINK, 1 , 1);
            zp = zeros(size(xp)) + Z_NKLINK;
        NECK_LINKER = [xp; yp; zp]';

            [xp, yp] = circus(R_tSHLINK, 1 , 1);
            zp = zeros(size(xp)) + Z_SHLINK;
        SHAFT_LINKER = [xp; yp; zp]';


    SPINE_ALL = [SPINE_PSD; SPINE_PERPSD; SPINE_SYNAPSE; SPINE_HEAD_UPPER; ... 
                 SPINE_HEAD_LOWER; SPINE_NECK; NECK_LINKER; SHAFT_LINKER];



  end


end



%% BUILD DENDRITE SKELLETON

function DENDRITE_ALL = get_dendrite(varargin)

%---------------- DENDRITIC SHAFT OBJECT -----------------%

    SCALE_FACTOR = 100;    % Mesh has issues below a certain scale

    DENDRITE_RADIUS = SCALE_FACTOR * 1.0;
    DENDRITE_LEFT = 0;

        [xp, zp] = circus(DENDRITE_RADIUS, 1 , 1);
        zp = zp - DENDRITE_RADIUS;
        yp = zeros(size(xp)) + DENDRITE_LEFT;
        DEN_L = [xp; yp; zp]';


    DEN_R = [];
    n = 50;    % number of rings
    nd = 10;    % distance between rings

    for nn = 1:n

        DENn = DEN_L; DENn(:,2) = DENn(:,2) + nn*nd;

        DEN_R = [DEN_R; DENn];

    end

    DENDRITE_ALL = [DEN_L; DEN_R];


end



%% HELPER FUNCTION TO CREATE CIRCLE POINTS

function [xp, yp, varargout] = circus(r,varargin)
% %------------------------------------------%
% % 3D Circle
% [xp yp zp] = circus(r,xc,yc,zc);
% figure
% [ph1] = plot3(xp, yp, zp);
% axis square;

ang=0:0.2:2*pi; 

if nargin == 1
	xc=0; yc=0;
	xp = r*cos(ang) + xc;
	yp = r*sin(ang) + yc;
elseif nargin == 2
	xc=varargin{1}; yc=0;
	xp = r*cos(ang) + xc;
	yp = r*sin(ang) + yc;
elseif nargin == 3
	xc=varargin{1}; yc=varargin{2};
	xp = r*cos(ang) + xc;
	yp = r*sin(ang) + yc;
elseif nargin == 4
	xc=varargin{1}; yc=varargin{2}; zc=varargin{3};
	xp = r*cos(ang) + xc;
	yp = r*sin(ang) + yc;
	zp = zeros(numel(xp)) + zc;
	varargout = {zp};
else
	disp('bad job');
end

end


end


