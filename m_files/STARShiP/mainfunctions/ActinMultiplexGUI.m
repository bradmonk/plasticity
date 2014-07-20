function varargout = ActinMultiplexGUI(varargin)
% ActinMultiplexGUI MATLAB code for ActinMultiplexGUI.fig
%      ActinMultiplexGUI, by itself, creates a new ActinMultiplexGUI or raises the existing
%      singleton*.
%
%      H = ActinMultiplexGUI returns the handle to a new ActinMultiplexGUI or the handle to
%      the existing singleton*.
%
%      ActinMultiplexGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ActinMultiplexGUI.M with the given input arguments.
%
%      ActinMultiplexGUI('Property','Value',...) creates a new ActinMultiplexGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ActinMultiplexGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ActinMultiplexGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ActinMultiplexGUI

% Last Modified by GUIDE v2.5 13-Feb-2014 18:46:44

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ActinMultiplexGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @ActinMultiplexGUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before ActinMultiplexGUI is made visible.
function ActinMultiplexGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ActinMultiplexGUI (see VARARGIN)

% Choose default command line output for ActinMultiplexGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% This sets up the initial plot - only do when we are invisible
% so window can get raised using ActinMultiplexGUI.
if strcmp(get(hObject,'Visible'),'off')
    % plot(rand(5));
end

% UIWAIT makes ActinMultiplexGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ActinMultiplexGUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

%=========================================================%
%=========================================================%
%=========================================================%
% --- Executes on button press in RunSAP.
function RunMultiplex_Callback(hObject, eventdata, handles)


%-------
Bon = str2double(get(handles.Bon,'String'));
Ron = str2double(get(handles.Ron,'String'));
Boff = str2double(get(handles.Boff,'String'));
Roff = str2double(get(handles.Roff,'String'));
Lon = str2double(get(handles.Lon,'String'));
Loff = str2double(get(handles.Loff,'String'));
%---
LBR = [Lon Loff Bon Boff Ron Roff];
%-------


%-------
Nsteps = str2double(get(handles.Nsteps,'String'));
dT = str2double(get(handles.dT,'String'));
DataR = str2double(get(handles.DataR,'String'));
ViewT = str2double(get(handles.ViewT,'String'));
PauseT = str2double(get(handles.PauseT,'String'));
LoopNum = str2double(get(handles.LoopNum,'String'));
%---
TIME = [Nsteps dT DataR ViewT PauseT LoopNum];
%-------


%-------
PSDsz = str2double(get(handles.PSDsz,'String'));
PSAsz = str2double(get(handles.PSAsz,'String'));
%---
SIZE = [PSDsz PSAsz];
%-------





%-------
doP1 = get(handles.doPlot,'Value');
doP2 = str2double(get(handles.PlotNum,'String'));
doP3 = get(handles.doFluorPlot,'Value');
doP4 = str2double(get(handles.FluorT,'String'));
%---
DOES = [doP1 doP2 doP3 doP4];
%-------



%-------
doRev = get(handles.doRev,'Value');
Rev = str2double(get(handles.Rev,'String'));
%---
REVA = [doRev Rev];
%-------





%-------
AMX1  = get(handles.loadActinTips,'Value');
AMX2  = get(handles.generateActinTips,'Value');
AMX3  = get(handles.ATfilename,'String');
AMX4  = str2double(get(handles.ActSteps,'String'));
AMX5  = str2double(get(handles.AMask,'String'));
AMX6  = str2double(get(handles.SMask,'String'));
AMX7  = str2double(get(handles.StartAct,'String'));
AMX8  = str2double(get(handles.ActUpdate,'String'));
AMX9  = str2double(get(handles.BranchAng,'String'));
AMX10  = str2double(get(handles.SPYneckXY,'String'));
AMX11  = str2double(get(handles.SPYheadZN,'String'));
AMX12  = str2double(get(handles.SPYheadZS,'String'));
AMX13  = str2double(get(handles.SPYheadX,'String'));
AMX14  = str2double(get(handles.SPYheadY,'String'));
AMX15  = str2double(get(handles.PSDproxy,'String'));
AMX16  = str2double(get(handles.molAct,'String'));
AMX17  = str2double(get(handles.ArpBR,'String'));
AMX18  = str2double(get(handles.CofR,'String'));
AMX19  = str2double(get(handles.SaveTipsAfter,'String'));
AMX20  = str2double(get(handles.SaveTipsRate,'String'));
AMX21  = str2double(get(handles.ARPsc,'String'));
AMX22  = str2double(get(handles.ARPmax,'String'));
AMX23  = str2double(get(handles.Asja,'String'));
AMX24  = str2double(get(handles.Asjb,'String'));
AMX25  = str2double(get(handles.Asjc,'String'));
AMX26  = str2double(get(handles.Asjd,'String'));
AMX27  = str2double(get(handles.PSDproxyXY,'String'));
AMX28  = get(handles.doClustering,'Value');
AMX29  = str2double(get(handles.ACTadd2new,'String'));
AMX30  = str2double(get(handles.ACTdelCofN,'String'));
AMX31  = str2double(get(handles.ArpON,'String'));
AMX32  = str2double(get(handles.ArpOFF,'String'));
AMX33  = str2double(get(handles.GArpN,'String'));
AMX34  = str2double(get(handles.CofScalar,'String'));
AMX35  = get(handles.doActCounts,'Value');
AMX36  = get(handles.doActLTP,'Value');
AMX37  = str2double(get(handles.AcLTPs,'String'));
AMX38  = str2double(get(handles.AcLTPe,'String'));
AMX39  = get(handles.DelOrigFils,'Value');
AMX40  = str2double(get(handles.DelOrigFilsT,'String'));
AMX41  = str2double(get(handles.startmonos,'String'));
AMX42  = str2double(get(handles.startfils,'String'));
AMX43  = str2double(get(handles.fZo,'String'));
AMX44  = str2double(get(handles.fZa,'String'));
AMX45  = str2double(get(handles.fXYo,'String'));
AMX46  = str2double(get(handles.fXYa,'String'));
AMX47  = str2double(get(handles.AdT,'String'));
AMX48  = str2double(get(handles.Aska,'String'));
AMX49  = str2double(get(handles.Askb,'String'));
AMX50  = str2double(get(handles.Askc,'String'));
AMX51  = str2double(get(handles.Askd,'String'));
AMX52  = str2double(get(handles.delOrate,'String'));
%---
AMX = {AMX1, AMX2, AMX3, AMX4, AMX5, AMX6, AMX7, AMX8, AMX9, AMX10,...
	   AMX11, AMX12, AMX13, AMX14, AMX15, AMX16, AMX17, AMX18, AMX19,...
	   AMX20, AMX21, AMX22, AMX23, AMX24, AMX25, AMX26, AMX27, AMX28,...
	   AMX29, AMX30, AMX31, AMX32, AMX33, AMX34, AMX35, AMX36, AMX37,...
	   AMX38, AMX39, AMX40, AMX41, AMX42, AMX43, AMX44, AMX45, AMX46,...
	   AMX47, AMX48, AMX49, AMX50, AMX51, AMX52};
%-------


%-------
AMS1  = str2double(get(handles.LivePlotM,'String'));
%--
AMS = {AMS1};
%-------


MSK1  = str2double(get(handles.GNpk,'String'));
MSK2  = str2double(get(handles.GNx0,'String'));
MSK3  = str2double(get(handles.GNy0,'String'));
MSK4  = str2double(get(handles.GNsd,'String'));
MSK5  = str2double(get(handles.GNnum,'String'));
MSK6  = str2double(get(handles.GNres,'String'));
%---
MSK = {MSK1, MSK2, MSK3, MSK4, MSK5, MSK6};
%-------






% handles.output = ActinMultiplex(LBR,TIME,SIZE,MODS,DOES,REVA,GLU,GT,GTab,doTs,AMX,MSK);
handles.output = ActinMultiplex(LBR,TIME,SIZE,DOES,REVA,AMX,MSK,AMS);



%=========================================================%
function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
function PrintMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to PrintMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
printdlg(handles.figure1)
function CloseMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
                     ['Close ' get(handles.figure1,'Name') '...'],...
                     'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end

delete(handles.figure1)
%=========================================================%






%----------Gaussian MSK  Callbacks--------------%
function GNpk(hObject, eventdata, handles)
GNpk = get(hObject, 'Value');
 
handles.guidata.GNpk = GNpk;
guidata(hObject,handles)
%---
function GNx0(hObject, eventdata, handles)
GNx0 = get(hObject, 'Value');
 
handles.guidata.GNx0 = GNx0;
guidata(hObject,handles)
%---
function GNy0(hObject, eventdata, handles)
GNy0 = get(hObject, 'Value');
 
handles.guidata.GNy0 = GNy0;
guidata(hObject,handles)
%---
function GNsd(hObject, eventdata, handles)
GNsd = get(hObject, 'Value');
 
handles.guidata.GNsd = GNsd;
guidata(hObject,handles)
%---
function GNnum(hObject, eventdata, handles)
GNnum = get(hObject, 'Value');
 
handles.guidata.GNnum = GNnum;
guidata(hObject,handles)
%---
function GNres(hObject, eventdata, handles)
GNres = get(hObject, 'Value');
 
handles.guidata.GNres = GNres;
guidata(hObject,handles)
%----------




%----------Actin Multiplex Callbacks--------------%
function loadActinTips(hObject, eventdata, handles)
loadActinTips = get(hObject, 'Value');
 
handles.guidata.loadActinTips = loadActinTips;
guidata(hObject,handles)
%---
function generateActinTips(hObject, eventdata, handles)
generateActinTips = get(hObject, 'Value');
 
handles.guidata.generateActinTips = generateActinTips;
guidata(hObject,handles)
%---
function ActSteps(hObject, eventdata, handles)
ActSteps = get(hObject, 'Value');
 
handles.guidata.ActSteps = ActSteps;
guidata(hObject,handles)
%---
function ATfilename(hObject, eventdata, handles)
ATfilename = get(hObject, 'Value');
 
handles.guidata.ATfilename = ATfilename;
guidata(hObject,handles)
%---
function AMask(hObject, eventdata, handles)
AMask = get(hObject, 'Value');
 
handles.guidata.AMask = AMask;
guidata(hObject,handles)
%---
function SMask(hObject, eventdata, handles)
SMask = get(hObject, 'Value');
 
handles.guidata.SMask = SMask;
guidata(hObject,handles)
%---
function StartAct(hObject, eventdata, handles)
StartAct = get(hObject, 'Value');
 
handles.guidata.StartAct = StartAct;
guidata(hObject,handles)
%---
function ActUpdate(hObject, eventdata, handles)
ActUpdate = get(hObject, 'Value');
 
handles.guidata.ActUpdate = ActUpdate;
guidata(hObject,handles)
%---
function BranchAng(hObject, eventdata, handles)
BranchAng = get(hObject, 'Value');
 
handles.guidata.BranchAng = BranchAng;
guidata(hObject,handles)
%---
function SPYheadX(hObject, eventdata, handles)
SPYheadX = get(hObject, 'Value');
 
handles.guidata.SPYheadX = SPYheadX;
guidata(hObject,handles)
%---
function SPYheadY(hObject, eventdata, handles)
SPYheadY = get(hObject, 'Value');
 
handles.guidata.SPYheadY = SPYheadY;
guidata(hObject,handles)
%---
function SPYneckXY(hObject, eventdata, handles)
SPYneckXY = get(hObject, 'Value');
 
handles.guidata.SPYneckXY = SPYneckXY;
guidata(hObject,handles)
%---
function SPYheadZN(hObject, eventdata, handles)
SPYheadZN = get(hObject, 'Value');
 
handles.guidata.SPYheadZN = SPYheadZN;
guidata(hObject,handles)
%---
function SPYheadZS(hObject, eventdata, handles)
SPYheadZS = get(hObject, 'Value');
 
handles.guidata.SPYheadZS = SPYheadZS;
guidata(hObject,handles)
%---
function PSDproxy(hObject, eventdata, handles)
PSDproxy = get(hObject, 'Value');
 
handles.guidata.PSDproxy = PSDproxy;
guidata(hObject,handles)
%---
function molAct(hObject, eventdata, handles)
molAct = get(hObject, 'Value');
 
handles.guidata.molAct = molAct;
guidata(hObject,handles)
%---
function ArpBR(hObject, eventdata, handles)
ArpBR = get(hObject, 'Value');
 
handles.guidata.ArpBR = ArpBR;
guidata(hObject,handles)
%---
function CofR(hObject, eventdata, handles)
CofR = get(hObject, 'Value');
 
handles.guidata.CofR = CofR;
guidata(hObject,handles)
%---
function SaveTipsAfter(hObject, eventdata, handles)
SaveTipsAfter = get(hObject, 'Value');
 
handles.guidata.SaveTipsAfter = SaveTipsAfter;
guidata(hObject,handles)
%---
function SaveTipsRate(hObject, eventdata, handles)
SaveTipsRate = get(hObject, 'Value');
 
handles.guidata.SaveTipsRate = SaveTipsRate;
guidata(hObject,handles)
%---
function ARPsc(hObject, eventdata, handles)
ARPsc = get(hObject, 'Value');
 
handles.guidata.ARPsc = ARPsc;
guidata(hObject,handles)
%---
function ARPmax(hObject, eventdata, handles)
ARPmax = get(hObject, 'Value');
 
handles.guidata.ARPmax = ARPmax;
guidata(hObject,handles)
%---
function Asja(hObject, eventdata, handles)
Asja = get(hObject, 'Value');
 
handles.guidata.Asja = Asja;
guidata(hObject,handles)
%---
function Asjb(hObject, eventdata, handles)
Asjb = get(hObject, 'Value');
 
handles.guidata.Asjb = Asjb;
guidata(hObject,handles)
%---
function Asjc(hObject, eventdata, handles)
Asjc = get(hObject, 'Value');
 
handles.guidata.Asjc = Asjc;
guidata(hObject,handles)
%---
function Asjd(hObject, eventdata, handles)
Asjd = get(hObject, 'Value');
 
handles.guidata.Asjd = Asjd;
guidata(hObject,handles)
%---
function PSDproxyXY(hObject, eventdata, handles)
PSDproxyXY = get(hObject, 'Value');
 
handles.guidata.PSDproxyXY = PSDproxyXY;
guidata(hObject,handles)
%---
function doClustering(hObject, eventdata, handles)
doClustering = get(hObject, 'Value');
 
handles.guidata.doClustering = doClustering;
guidata(hObject,handles)
%---
function ACTadd2new(hObject, eventdata, handles)
ACTadd2new = get(hObject, 'Value');
 
handles.guidata.ACTadd2new = ACTadd2new;
guidata(hObject,handles)
%---
function ACTdelCofN(hObject, eventdata, handles)
ACTdelCofN = get(hObject, 'Value');
 
handles.guidata.ACTdelCofN = ACTdelCofN;
guidata(hObject,handles)
%---
function ArpON(hObject, eventdata, handles)
ArpON = get(hObject, 'Value');
 
handles.guidata.ArpON = ArpON;
guidata(hObject,handles)
%---
function ArpOFF(hObject, eventdata, handles)
ArpOFF = get(hObject, 'Value');
 
handles.guidata.ArpOFF = ArpOFF;
guidata(hObject,handles)
%---
function GArpN(hObject, eventdata, handles)
GArpN = get(hObject, 'Value');
 
handles.guidata.GArpN = GArpN;
guidata(hObject,handles)
%---
function CofScalar(hObject, eventdata, handles)
CofScalar = get(hObject, 'Value');
 
handles.guidata.CofScalar = CofScalar;
guidata(hObject,handles)
%---
function doActCounts(hObject, eventdata, handles)
doActCounts = get(hObject, 'Value');
 
handles.guidata.doActCounts = doActCounts;
guidata(hObject,handles)
%---
function doActLTP(hObject, eventdata, handles)
doActLTP = get(hObject, 'Value');
 
handles.guidata.doActLTP = doActLTP;
guidata(hObject,handles)
%---
function AcLTPs(hObject, eventdata, handles)
AcLTPs = get(hObject, 'Value');
 
handles.guidata.AcLTPs = AcLTPs;
guidata(hObject,handles)
%---
function AcLTPe(hObject, eventdata, handles)
AcLTPe = get(hObject, 'Value');
 
handles.guidata.AcLTPe = AcLTPe;
guidata(hObject,handles)
%---
function DelOrigFils(hObject, eventdata, handles)
DelOrigFils = get(hObject, 'Value');
 
handles.guidata.DelOrigFils = DelOrigFils;
guidata(hObject,handles)
%---
function DelOrigFilsT(hObject, eventdata, handles)
DelOrigFilsT = get(hObject, 'Value');
 
handles.guidata.DelOrigFilsT = DelOrigFilsT;
guidata(hObject,handles)
%---
function startmonos(hObject, eventdata, handles)
startmonos = get(hObject, 'Value');
 
handles.guidata.startmonos = startmonos;
guidata(hObject,handles)
%---
function startfils(hObject, eventdata, handles)
startfils = get(hObject, 'Value');
 
handles.guidata.startfils = startfils;
guidata(hObject,handles)
%---
function fZo(hObject, eventdata, handles)
fZo = get(hObject, 'Value');
 
handles.guidata.fZo = fZo;
guidata(hObject,handles)
%---
function fZa(hObject, eventdata, handles)
fZa = get(hObject, 'Value');
 
handles.guidata.fZa = fZa;
guidata(hObject,handles)
%---
function fXYo(hObject, eventdata, handles)
fXYo = get(hObject, 'Value');
 
handles.guidata.fXYo = fXYo;
guidata(hObject,handles)
%---
function fXYa(hObject, eventdata, handles)
fXYa = get(hObject, 'Value');
 
handles.guidata.fXYa = fXYa;
guidata(hObject,handles)
%---
function AdT(hObject, eventdata, handles)
AdT = get(hObject, 'Value');
 
handles.guidata.AdT = AdT;
guidata(hObject,handles)
%---
function Aska(hObject, eventdata, handles)
Aska = get(hObject, 'Value');
 
handles.guidata.Aska = Aska;
guidata(hObject,handles)
%---
function Askb(hObject, eventdata, handles)
Askb = get(hObject, 'Value');
 
handles.guidata.Askb = Askb;
guidata(hObject,handles)
%---
function Askc(hObject, eventdata, handles)
Askc = get(hObject, 'Value');
 
handles.guidata.Askc = Askc;
guidata(hObject,handles)
%---
function Askd(hObject, eventdata, handles)
Askd = get(hObject, 'Value');
 
handles.guidata.Askd = Askd;
guidata(hObject,handles)
%---
function delOrate(hObject, eventdata, handles)
delOrate = get(hObject, 'Value');
 
handles.guidata.delOrate = delOrate;
guidata(hObject,handles)
%---
function LivePlotM(hObject, eventdata, handles)
LivePlotM = get(hObject, 'Value');
 
handles.guidata.LivePlotM = LivePlotM;
guidata(hObject,handles)
%---




%{

LivePlotM

function XXXXX(hObject, eventdata, handles)
XXXXX = get(hObject, 'Value');
 
handles.guidata.XXXXX = XXXXX;
guidata(hObject,handles)
%---



function XXXXX(hObject, eventdata, handles)
XXXXX = get(hObject, 'Value');
 
handles.guidata.XXXXX = XXXXX;
guidata(hObject,handles)
%---




function XXXXX(hObject, eventdata, handles)
XXXXX = get(hObject, 'Value');
 
handles.guidata.XXXXX = XXXXX;
guidata(hObject,handles)
%---



function XXXXX(hObject, eventdata, handles)
XXXXX = get(hObject, 'Value');
 
handles.guidata.XXXXX = XXXXX;
guidata(hObject,handles)
%---

%}

















%----------Time Value Callbacks--------------%
%---
function do1h(hObject, eventdata, handles)
do1h = get(hObject, 'Value');
 
handles.guidata.do1h = do1h;
guidata(hObject,handles)
%---
function do4h(hObject, eventdata, handles)
do4h = get(hObject, 'Value');
 
handles.guidata.do4h = do4h;
guidata(hObject,handles)
%---
function do1d(hObject, eventdata, handles)
do1d = get(hObject, 'Value');
 
handles.guidata.do1d = do1d;
guidata(hObject,handles)
%---
function do1m(hObject, eventdata, handles)
do1m = get(hObject, 'Value');
 
handles.guidata.do1m = do1m;
guidata(hObject,handles)
%---
function do1y(hObject, eventdata, handles)
do1y = get(hObject, 'Value');
 
handles.guidata.do1y = do1y;
guidata(hObject,handles)
%--------
function do1hT(hObject, eventdata, handles)
do1hT = get(hObject, 'Value');
 
handles.guidata.do1hT = do1hT;
guidata(hObject,handles)
%---
function do4hT(hObject, eventdata, handles)
do4hT = get(hObject, 'Value');
 
handles.guidata.do4hT = do4hT;
guidata(hObject,handles)
%---
function do1dT(hObject, eventdata, handles)
do1dT = get(hObject, 'Value');
 
handles.guidata.do1dT = do1dT;
guidata(hObject,handles)
%---
function do1mT(hObject, eventdata, handles)
do1mT = get(hObject, 'Value');
 
handles.guidata.do1mT = do1mT;
guidata(hObject,handles)
%---
function do1yT(hObject, eventdata, handles)
do1yT = get(hObject, 'Value');
 
handles.guidata.do1yT = do1yT;
guidata(hObject,handles)


%----------Numeric Value Callbacks--------------%
function GTon(hObject, eventdata, handles)
GTon = str2double(get(hObject, 'String'));
if isnan(GTon)
    GTon = 10;
end
handles.guidata.GTon = GTon;
guidata(hObject,handles)
%---
function GToff(hObject, eventdata, handles)
GToff = str2double(get(hObject, 'String'));
if isnan(GToff)
    GT1off = -5;
end
handles.guidata.GToff = GToff;
guidata(hObject,handles)
%---
function LTPv(hObject, eventdata, handles)
LTPv = str2double(get(hObject, 'String'));
if isnan(LTPv)
    LTPv = .001;
end
 
% Save the new psd1D value
handles.guidata.LTPv = LTPv;
guidata(hObject,handles)
%---
function GMASK(hObject, eventdata, handles)
GMASK = get(hObject, 'Data');
if isnan(GMASK)
    GMASK = logical([0 1 0; 1 1 1; 0 1 0]);
end
handles.guidata.GMASK = GMASK;
guidata(hObject,handles)
%---
function gM(hObject, eventdata, handles)
gM = get(hObject, 'Value');
 
handles.guidata.gM = gM;
guidata(hObject,handles)

%----------Numeric Value Callbacks--------------%
function Bon(hObject, eventdata, handles)
Bon = str2double(get(hObject, 'String'));
if isnan(Bon)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

handles.guidata.Bon = Bon;
guidata(hObject,handles)
%---
function Ron(hObject, eventdata, handles)
Ron = str2double(get(hObject, 'String'));
if isnan(Ron)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

handles.guidata.Ron = Ron;
guidata(hObject,handles)
%---
function Boff(hObject, eventdata, handles)
Boff = str2double(get(hObject, 'String'));
if isnan(Boff)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

handles.guidata.Boff = Boff;
guidata(hObject,handles)
%---
function Roff(hObject, eventdata, handles)
Roff = str2double(get(hObject, 'String'));
if isnan(Roff)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

handles.guidata.Roff = Roff;
guidata(hObject,handles)
%---
function Lon(hObject, eventdata, handles)
Lon = str2double(get(hObject, 'String'));
if isnan(Bon)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

handles.guidata.Lon = Lon;
guidata(hObject,handles)
%---
function Loff(hObject, eventdata, handles)
Loff = str2double(get(hObject, 'String'));
if isnan(Ron)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

handles.guidata.Loff = Loff;
guidata(hObject,handles)
%---
function Nsteps(hObject, eventdata, handles)
Nsteps = str2double(get(hObject, 'String'));
if isnan(Nsteps)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

handles.guidata.Nsteps = Nsteps;
guidata(hObject,handles)
%---
function dT(hObject, eventdata, handles)
dT = str2double(get(hObject, 'String'));
if isnan(dT)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

handles.guidata.dT = dT;
guidata(hObject,handles)
%---
function G1RT(hObject, eventdata, handles)
G1RT = str2double(get(hObject, 'String'));
if isnan(G1RT)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

handles.guidata.G1RT = G1RT;
guidata(hObject,handles)
%---
function G1ST(hObject, eventdata, handles)
G1ST = str2double(get(hObject, 'String'));
if isnan(G1ST)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

handles.guidata.G1ST = G1ST;
guidata(hObject,handles)
%---
function AMPARrate(hObject, eventdata, handles)
AMPARrate = str2double(get(hObject, 'String'));
if isnan(AMPARrate)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

handles.guidata.AMPARrate = AMPARrate;
guidata(hObject,handles)
%---
function PSDsz(hObject, eventdata, handles)
PSDsz = str2double(get(hObject, 'String'));
if isnan(PSDsz)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

handles.guidata.PSDsz = PSDsz;
guidata(hObject,handles)
%---
function PSAsz(hObject, eventdata, handles)
PSAsz = str2double(get(hObject, 'String'));
if isnan(PSAsz)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

handles.guidata.PSAsz = PSAsz;
guidata(hObject,handles)
%---
function DataR(hObject, eventdata, handles)
DataR = str2double(get(hObject, 'String'));
if isnan(DataR)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

handles.guidata.DataR = DataR;
guidata(hObject,handles)
%---
function ViewT(hObject, eventdata, handles)
ViewT = str2double(get(hObject, 'String'));
if isnan(ViewT)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

handles.guidata.ViewT = ViewT;
guidata(hObject,handles)
%---
function PauseT(hObject, eventdata, handles)
PauseT = str2double(get(hObject, 'String'));
if isnan(PauseT)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

handles.guidata.PauseT = PauseT;
guidata(hObject,handles)
%---
function FluorT(hObject, eventdata, handles)
FluorT = str2double(get(hObject, 'String'));
if isnan(FluorT)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

handles.guidata.FluorT = FluorT;
guidata(hObject,handles)
%---
function PlotNum(hObject, eventdata, handles)
PlotNum = str2double(get(hObject, 'String'));
if isnan(PlotNum)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

handles.guidata.PlotNum = PlotNum;
guidata(hObject,handles)
%---
function Rev(hObject, eventdata, handles)
Rev = str2double(get(hObject, 'String'));
if isnan(Rev)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

handles.guidata.Rev = Rev;
guidata(hObject,handles)
%---
function LoopNum(hObject, eventdata, handles)
LoopNum = str2double(get(hObject, 'String'));
if isnan(LoopNum)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

handles.guidata.LoopNum = LoopNum;
guidata(hObject,handles)
%---
function AMPARnum(hObject, eventdata, handles)
AMPARnum = str2double(get(hObject, 'String'));
if isnan(AMPARnum)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

handles.guidata.AMPARnum = AMPARnum;
guidata(hObject,handles)
%---
function ORM(hObject, eventdata, handles)
ORM = get(hObject, 'Value');
 
handles.guidata.ORM = ORM;
guidata(hObject,handles)
%---
function OffRM(hObject, eventdata, handles)
OffRM = get(hObject, 'Value');
 
handles.guidata.OffRM = OffRM;
guidata(hObject,handles)
%---
function OnRM(hObject, eventdata, handles)
OnRM = get(hObject, 'Value');
 
handles.guidata.OnRM = OnRM;
guidata(hObject,handles)
%---
function doAMPAR(hObject, eventdata, handles)
doAMPAR = get(hObject, 'Value');
 
handles.guidata.doAMPAR = doAMPAR;
guidata(hObject,handles)
%---
function doPlot(hObject, eventdata, handles)
doPlot = get(hObject, 'Value');
 
handles.guidata.doPlot = doPlot;
guidata(hObject,handles)
%---
function doRev(hObject, eventdata, handles)
doRev = get(hObject, 'Value');
 
handles.guidata.doRev = doRev;
guidata(hObject,handles)
%---
function doFluorPlot(hObject, eventdata, handles)
doRev = get(hObject, 'Value');
 
handles.guidata.doRev = doRev;
guidata(hObject,handles)
%---

 
%----------Logical Callbacks--------------%
%{
function doSpike(hObject, eventdata, handles)
doSpike = get(hObject, 'Value');
 
handles.guidata.doSpike = doSpike;
guidata(hObject,handles)
%}
 
 
 
 
 
 
 
%=========================================================%
%=========================================================%
% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
 
% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
     set(hObject,'BackgroundColor','white');
end
 
set(hObject, 'String', {'plot(rand(5))', 'plot(sin(1:0.01:25))', 'bar(1:.5:10)', 'plot(membrane)', 'surf(peaks)'});

function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 
% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double

% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
 
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 
% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double
 
 
% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
 
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
 
function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 
% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double

% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
 
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
 
 
 
function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 
% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double
 
 
% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
 
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
 
%=========================================================%
%=========================================================%








