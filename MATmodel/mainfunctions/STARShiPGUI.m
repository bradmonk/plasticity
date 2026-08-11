function varargout = STARShiPGUI(varargin)
% BradsModelGUI MATLAB code for BradsModelGUI.fig
%      BradsModelGUI, by itself, creates a new BradsModelGUI or raises the existing
%      singleton*.
%
%      H = BradsModelGUI returns the handle to a new BradsModelGUI or the handle to
%      the existing singleton*.
%
%      BradsModelGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BradsModelGUI.M with the given input arguments.
%
%      BradsModelGUI('Property','Value',...) creates a new BradsModelGUI or raises
%      the existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before BradsModelGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to BradsModelGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help BradsModelGUI

% Last Modified by GUIDE v2.5 06-Oct-2013 17:23:26

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @BradsModelGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @BradsModelGUI_OutputFcn, ...
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

% --- Executes just before BradsModelGUI is made visible.
function BradsModelGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to BradsModelGUI (see VARARGIN)

% Choose default command line output for BradsModelGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

initialize_gui(hObject, handles, false);

% UIWAIT makes BradsModelGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = BradsModelGUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%  reset.
function reset_Callback(hObject, eventdata, handles)
% hObject    handle to reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

initialize_gui(gcbf, handles, true);


function initialize_gui(fig_handle, handles, isreset)

% set(handles.GluR2dots, 'String', handles.GluR2dots);
% set(handles.GluR1dots,  'String', handles.GluR1dots);
% set(handles.GraphTime,  'String', handles.GraphTime);
% set(handles.AllowedTime,  'String', handles.AllowedTime);
% set(handles.spiD,  'String', handles.spiD);
% set(handles.psd2D,  'String', handles.psd2D);
% set(handles.psd1D,  'String', handles.psd1D);

% Update handles structure
guidata(handles.figure1, handles);

%========================================================================%
%------------------------------------------------------------------------%
%========================================================================%

%  RUN SIMULATION STARShiP
function RunSimulation_Callback(hObject, eventdata, handles)



%-------
dot1 = str2double(get(handles.GluR1dots,'String'));
dot2 = str2double(get(handles.GluR2dots,'String'));
dot3 = str2double(get(handles.Steps,'String'));
dot4 = str2double(get(handles.TimeStep,'String'));
dot5 = str2double(get(handles.Scale,'String'));
dot6 = str2double(get(handles.doLoops,'String'));
dot7 = str2double(get(handles.MplotUp,'String'));
dot8 = str2double(get(handles.SplotUp,'String'));
%---
dot = [dot1 dot2 dot3 dot4 dot5 dot6 dot7 dot8];
%-------



%-------
dr1 = str2double(get(handles.esDGR1,'String'));
dr2 = str2double(get(handles.spineDGR1,'String'));
dr3 = str2double(get(handles.PERIDGR1,'String'));
dr4 = str2double(get(handles.PSDDGR1,'String'));
dr5 = str2double(get(handles.esDGR2,'String'));
dr6 = str2double(get(handles.spineDGR2,'String'));
dr7 = str2double(get(handles.PERIDGR2,'String'));
dr8 = str2double(get(handles.PSDDGR2,'String'));
%---
dr = [dr1 dr2 dr3 dr4 dr5 dr6 dr7 dr8];
%-------



%-------
um1 = str2double(get(handles.denWidthX,'String'));
um2 = str2double(get(handles.denHeightY,'String'));
um3 = str2double(get(handles.PSD1um,'String'));
um4 = str2double(get(handles.PSD2um,'String'));
um5 = str2double(get(handles.PERI1um,'String'));
um6 = str2double(get(handles.PERI2um,'String'));
%---
um = [um1 um2 um3 um4 um5 um6];
%-------



%-------
sap11 = str2double(get(handles.SAPdTPSD1,'String'));
sap12 = str2double(get(handles.SAPdTPSD2,'String'));
sap19 = get(handles.doLTPS1,'Value');
sap20 = get(handles.doLTPS2,'Value');
sap21 = str2double(get(handles.LonS1,'String'));
sap22 = str2double(get(handles.BonS1,'String'));
sap23 = str2double(get(handles.RonS1,'String'));
sap24 = str2double(get(handles.LoffS1,'String'));
sap25 = str2double(get(handles.BoffS1,'String'));
sap26 = str2double(get(handles.RoffS1,'String'));
sap27 = str2double(get(handles.LonS2,'String'));
sap28 = str2double(get(handles.BonS2,'String'));
sap29 = str2double(get(handles.RonS2,'String'));
sap30 = str2double(get(handles.LoffS2,'String'));
sap31 = str2double(get(handles.BoffS2,'String'));
sap32 = str2double(get(handles.RoffS2,'String'));
%---
sap = [0 0 0 0 0 0 0 0 0 0 sap11 sap12 0 0 0 0 0 0 sap19 sap20...
	   sap21 sap22 sap23 sap24 sap25 sap26 sap27 sap28 sap29 sap30 sap31 sap32];
%-------


%-------
atn1 = get(handles.AS1doup,'Value');
atn2 = get(handles.AS2doup,'Value');
atn3 = str2double(get(handles.AS1updt,'String'));
atn4 = str2double(get(handles.AS2updt,'String'));
atn5 = get(handles.FileATdataS1,'String');
atn6 = get(handles.FileATdataS2,'String');
atn7 = str2double(get(handles.ActS1scell,'String'));
atn8 = str2double(get(handles.ActS2scell,'String'));
%---
atn = {atn1,atn2,atn3,atn4,atn5,atn6,atn7,atn8};
%-------


%-------
doUse1 = get(handles.useGluR1,'Value');
doUse2 = get(handles.useGluR2,'Value');
doUse3 = get(handles.runSAPPSD1,'Value');
doUse4 = get(handles.runSAPPSD2,'Value');
%---
doUse = [doUse1 doUse2 doUse3 doUse4];
%-------



%-------
doRun1 = get(handles.run2Dplot,'Value');
doRun2 = get(handles.run3Dplot,'Value');
doRun3 = get(handles.runMSDtest,'Value');
doRun4 = 0;
doRun5 = get(handles.runMSDpopup,'Value');		% POPUP
doRun6 = 0;
doRun8 = get(handles.doFieldFig,'Value');
doRun9 = get(handles.doSlotColormap,'Value');
doRun10 = get(handles.doSpike,'Value');
doRun12 = get(handles.KBrd,'Value');
%---
doRun = [doRun1 doRun2 doRun3 doRun4 doRun5 doRun6 1 doRun8 doRun9 doRun10 0 doRun12];
%-------





handles.output = STARShiP(dot,dr,um,sap,doUse,doRun,atn);


% handles.output = ReDiClus(dot,dr,um,sap,hr,ko,doUse,doRun,doKo,box,slt,stky);
%========================================================================%
%------------------------------------------------------------------------%
%========================================================================%


%=========================================================%
function GluR2dots(hObject, eventdata, handles)

GluR2dots = str2double(get(hObject, 'String'));
if isnan(GluR2dots)
    set(hObject, 'String', 80);
    errordlg('Input must be a number','Error');
end

% Save the new GluR2dots value
handles.guidata.GluR2dots = GluR2dots;
guidata(hObject,handles)
%---
function GraphTime(hObject, eventdata, handles)
GraphTime = str2double(get(hObject, 'String'));
if isnan(GraphTime)
    set(hObject, 'String', 30);
    errordlg('Input must be a number','Error');
end

% Save the new GraphTime value
handles.guidata.GraphTime = GraphTime;
guidata(hObject,handles)
%---
function AllowedTime(hObject, eventdata, handles)
AllowedTime = str2double(get(hObject, 'String'));
if isnan(AllowedTime)
    set(hObject, 'String', 200);
    errordlg('Input must be a number','Error');
end

% Save the new AllowedTime value
handles.guidata.AllowedTime = AllowedTime;
guidata(hObject,handles)
%---
function GluR1dots(hObject, eventdata, handles)
GluR1dots = str2double(get(hObject, 'String'));
if isnan(GluR1dots)
    set(hObject, 'String', 20);
    errordlg('Input must be a number','Error');
end

% Save the new GluR1dots value
handles.guidata.GluR1dots = GluR1dots;
guidata(hObject,handles)
%---
function spiD(hObject, eventdata, handles)
spiD = str2double(get(hObject, 'String'));
if isnan(spiD)
    set(hObject, 'String', 0.1);
    errordlg('Input must be a number','Error');
end

% Save the new spiD value
handles.guidata.spiD = spiD;
guidata(hObject,handles)
%---
function psd2D(hObject, eventdata, handles)
psd2D = str2double(get(hObject, 'String'));
if isnan(psd2D)
    set(hObject, 'String', 0.01);
    errordlg('Input must be a number','Error');
end

% Save the new psd2D value
handles.guidata.psd2D = psd2D;
guidata(hObject,handles)
%---
function psd1D(hObject, eventdata, handles)
psd1D = str2double(get(hObject, 'String'));
if isnan(psd1D)
    set(hObject, 'String', 0.001);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.psd1D = psd1D;
guidata(hObject,handles)
%---
function NumLoops(hObject, eventdata, handles)
NumLoops = str2double(get(hObject, 'String'));
if isnan(NumLoops)
    set(hObject, 'String', 5);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.NumLoops = NumLoops;
guidata(hObject,handles)
%---
function SteadyState(hObject, eventdata, handles)
% Hint: get(hObject,'Value') returns toggle state of SteadyState

SteadyState = get(hObject, 'Value');

% Save the new SteadyState value
handles.guidata.SteadyState = SteadyState;
guidata(hObject,handles)
%---
function doLiveSim(hObject, eventdata, handles)
% Hint: get(hObject,'Value') returns toggle state of doLiveSim

doLiveSim = get(hObject, 'Value');

% Save the new doLiveSim value
handles.guidata.doLiveSim = doLiveSim;
guidata(hObject,handles)
%=========================================================%
function doLoops(hObject, eventdata, handles)
doLoops = str2double(get(hObject, 'String'));
if isnan(doLoops)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.doLoops = doLoops;
guidata(hObject,handles)
%---
function Steps(hObject, eventdata, handles)
Steps = str2double(get(hObject, 'String'));
if isnan(Steps)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.Steps = Steps;
guidata(hObject,handles)
%---
function TimeStep(hObject, eventdata, handles)
TimeStep = str2double(get(hObject, 'String'));
if isnan(TimeStep)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.TimeStep = TimeStep;
guidata(hObject,handles)
%---
function Scale(hObject, eventdata, handles)
Scale = str2double(get(hObject, 'String'));
if isnan(Scale)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.Scale = Scale;
guidata(hObject,handles)
%=========================================================%
function runMSDtest(hObject, eventdata, handles)
runMSDtest = get(hObject, 'Value');

% Save the new SteadyState value
handles.guidata.runMSDtest = runMSDtest;
guidata(hObject,handles)
%---
function runTraceSingleDot(hObject, eventdata, handles)
runTraceSingleDot = get(hObject, 'Value');

% Save the new SteadyState value
handles.guidata.runTraceSingleDot = runTraceSingleDot;
guidata(hObject,handles)
%---
function runMSDpopup(hObject, eventdata, handles)
runMSDpopup = get(hObject, 'Value');

% Save the new SteadyState value
handles.guidata.runMSDpopup = runMSDpopup;
guidata(hObject,handles)
%---
function run2Dplot(hObject, eventdata, handles)
run2Dplot = get(hObject, 'Value');

% Save the new SteadyState value
handles.guidata.run2Dplot = run2Dplot;
guidata(hObject,handles)
%---
function run3Dplot(hObject, eventdata, handles)
run3Dplot = get(hObject, 'Value');

% Save the new SteadyState value
handles.guidata.run3Dplot = run3Dplot;
guidata(hObject,handles)
%=========================================================%
function SAPdotsPSD1(hObject, eventdata, handles)
SAPdotsPSD1 = str2double(get(hObject, 'String'));
if isnan(SAPdotsPSD1)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.SAPdotsPSD1 = SAPdotsPSD1;
guidata(hObject,handles)
%---
function SAPdotsPSD2(hObject, eventdata, handles)
SAPdotsPSD2 = str2double(get(hObject, 'String'));
if isnan(SAPdotsPSD2)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.SAPdotsPSD2 = SAPdotsPSD2;
guidata(hObject,handles)
%---
function runSAPPSD1(hObject, eventdata, handles)
runSAPPSD1 = get(hObject, 'Value');

% Save the new SteadyState value
handles.guidata.runSAPPSD1 = runSAPPSD1;
guidata(hObject,handles)
%---
function runSAPPSD2(hObject, eventdata, handles)
runSAPPSD2 = get(hObject, 'Value');

% Save the new SteadyState value
handles.guidata.runSAPPSD2 = runSAPPSD2;
guidata(hObject,handles)
%---
function doDynamicLeP1(hObject, eventdata, handles)
doDynamicLeP1 = get(hObject, 'Value');

% Save the new SteadyState value
handles.guidata.doDynamicLeP1 = doDynamicLeP1;
guidata(hObject,handles)
%---
function doDynamicLeP2(hObject, eventdata, handles)
doDynamicLeP2 = get(hObject, 'Value');

% Save the new SteadyState value
handles.guidata.doDynamicLeP2 = doDynamicLeP2;
guidata(hObject,handles)
%---
function doLTPS1(hObject, eventdata, handles)
doLTPS1 = get(hObject, 'Value');

% Save the new SteadyState value
handles.guidata.doLTPS1 = doLTPS1;
guidata(hObject,handles)
%---
function doLTPS2(hObject, eventdata, handles)
doLTPS2 = get(hObject, 'Value');

% Save the new SteadyState value
handles.guidata.doLTPS2 = doLTPS2;
guidata(hObject,handles)
%---
function denWidthX(hObject, eventdata, handles)
denWidthX = str2double(get(hObject, 'String'));
if isnan(denWidthX)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.denWidthX = denWidthX;
guidata(hObject,handles)
%---
function denHeightY(hObject, eventdata, handles)
denHeightY = str2double(get(hObject, 'String'));
if isnan(denHeightY)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.denHeightY = denHeightY;
guidata(hObject,handles)
%---
function PSD1um(hObject, eventdata, handles)
PSD1um = str2double(get(hObject, 'String'));
if isnan(PSD1um)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.PSD1um = PSD1um;
guidata(hObject,handles)
%---
function PSD2um(hObject, eventdata, handles)
PSD2um = str2double(get(hObject, 'String'));
if isnan(PSD2um)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.PSD2um = PSD2um;
guidata(hObject,handles)
%---
function PERI1um(hObject, eventdata, handles)
PERI1um = str2double(get(hObject, 'String'));
if isnan(PERI1um)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.PERI1um = PERI1um;
guidata(hObject,handles)
%---
function PERI2um(hObject, eventdata, handles)
PERI2um = str2double(get(hObject, 'String'));
if isnan(PERI2um)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.PERI2um = PERI2um;
guidata(hObject,handles)
%---
function spineDGR1(hObject, eventdata, handles)
spineDGR1 = str2double(get(hObject, 'String'));
if isnan(spineDGR1)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.spineDGR1 = spineDGR1;
guidata(hObject,handles)
%---
function PERIDGR1(hObject, eventdata, handles)
PERIDGR1 = str2double(get(hObject, 'String'));
if isnan(PERIDGR1)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.PERIDGR1 = PERIDGR1;
guidata(hObject,handles)
%---
function PSDDGR1(hObject, eventdata, handles)
PSDDGR1 = str2double(get(hObject, 'String'));
if isnan(PSDDGR1)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.PSDDGR1 = PSDDGR1;
guidata(hObject,handles)
%---
function esDGR1(hObject, eventdata, handles)
esDGR1 = str2double(get(hObject, 'String'));
if isnan(esDGR1)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.esDGR1 = esDGR1;
guidata(hObject,handles)
%---
function spineDGR2(hObject, eventdata, handles)
spineDGR2 = str2double(get(hObject, 'String'));
if isnan(spineDGR2)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.spineDGR2 = spineDGR2;
guidata(hObject,handles)
%---
function PERIDGR2(hObject, eventdata, handles)
PERIDGR2 = str2double(get(hObject, 'String'));
if isnan(PERIDGR2)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.PERIDGR2 = PERIDGR2;
guidata(hObject,handles)
%---
function PSDDGR2(hObject, eventdata, handles)
PSDDGR2 = str2double(get(hObject, 'String'));
if isnan(PSDDGR2)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.PSDDGR2 = PSDDGR2;
guidata(hObject,handles)
%---
function esDGR2(hObject, eventdata, handles)
esDGR2 = str2double(get(hObject, 'String'));
if isnan(esDGR2)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.esDGR2 = esDGR2;
guidata(hObject,handles)
%---
function useGluR1(hObject, eventdata, handles)
useGluR1 = get(hObject, 'Value');

% Save the new SteadyState value
handles.guidata.useGluR1 = useGluR1;
guidata(hObject,handles)
%---
function useGluR2(hObject, eventdata, handles)
useGluR2 = get(hObject, 'Value');

% Save the new SteadyState value
handles.guidata.useGluR2 = useGluR2;
guidata(hObject,handles)
%=========================================================%
function doFieldFig(hObject, eventdata, handles)
doFieldFig = get(hObject, 'Value');

handles.guidata.doFieldFig = doFieldFig;
guidata(hObject,handles)
%---
function doSlotColormap(hObject, eventdata, handles)
doSlotColormap = get(hObject, 'Value');

handles.guidata.doSlotColormap = doSlotColormap;
guidata(hObject,handles)
%=========================================================%
function LonS1(hObject, eventdata, handles)
LonS1 = str2double(get(hObject, 'String'));
if isnan(LonS1)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

handles.guidata.LonS1 = LonS1;
guidata(hObject,handles)
%---
function BonS1(hObject, eventdata, handles)
BonS1 = get(hObject, 'Value');

handles.guidata.BonS1 = BonS1;
guidata(hObject,handles)
%---
function RonS1(hObject, eventdata, handles)
RonS1 = str2double(get(hObject, 'String'));
if isnan(RonS1)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

handles.guidata.RonS1 = RonS1;
guidata(hObject,handles)
%---
function LoffS1(hObject, eventdata, handles)
LoffS1 = str2double(get(hObject, 'String'));
if isnan(LoffS1)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

handles.guidata.LoffS1 = LoffS1;
guidata(hObject,handles)
%---
function BoffS1(hObject, eventdata, handles)
BoffS1 = get(hObject, 'Value');

handles.guidata.BoffS1 = BoffS1;
guidata(hObject,handles)
%---
function RoffS1(hObject, eventdata, handles)
RoffS1 = str2double(get(hObject, 'String'));
if isnan(RoffS1)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

handles.guidata.RoffS1 = RoffS1;
guidata(hObject,handles)
%---
function LonS2(hObject, eventdata, handles)
LonS2 = str2double(get(hObject, 'String'));
if isnan(LonS2)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end
 
handles.guidata.LonS2 = LonS2;
guidata(hObject,handles)
%---
function BonS2(hObject, eventdata, handles)
BonS2 = get(hObject, 'Value');

handles.guidata.BonS2 = BonS2;
guidata(hObject,handles)
%---
function RonS2(hObject, eventdata, handles)
RonS2 = str2double(get(hObject, 'String'));
if isnan(RonS2)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end
 
handles.guidata.RonS2 = RonS2;
guidata(hObject,handles)
%---
function LoffS2(hObject, eventdata, handles)
LoffS2 = str2double(get(hObject, 'String'));
if isnan(LoffS2)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end
 
handles.guidata.LoffS2 = LoffS2;
guidata(hObject,handles)
%---
function BoffS2(hObject, eventdata, handles)
BoffS2 = get(hObject, 'Value');

handles.guidata.BoffS2 = BoffS2;
guidata(hObject,handles)
%---
function RoffS2(hObject, eventdata, handles)
RoffS2 = str2double(get(hObject, 'String'));
if isnan(RoffS2)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end
 
handles.guidata.RoffS2 = RoffS2;
guidata(hObject,handles)
%---
function SAPdTPSD1(hObject, eventdata, handles)
SAPdTPSD1 = str2double(get(hObject, 'String'));
if isnan(SAPdTPSD1)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end
 
% Save the new psd1D value
handles.guidata.SAPdTPSD1 = SAPdTPSD1;
guidata(hObject,handles)
%---
function SAPdTPSD2(hObject, eventdata, handles)
SAPdTPSD2 = str2double(get(hObject, 'String'));
if isnan(SAPdTPSD2)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end
 
% Save the new psd1D value
handles.guidata.SAPdTPSD2 = SAPdTPSD2;
guidata(hObject,handles)
%=========================================================%
function doSpike(hObject, eventdata, handles)
doSpike = get(hObject, 'Value');

handles.guidata.doSpike = doSpike;
guidata(hObject,handles)
%---
function doProfile(hObject, eventdata, handles)
doProfile = get(hObject, 'Value');

handles.guidata.doProfile = doProfile;
guidata(hObject,handles)
%=========================================================%
function AS1doup(hObject, eventdata, handles)
AS1doup = get(hObject, 'Value');

handles.guidata.AS1doup = AS1doup;
guidata(hObject,handles)
%---
function AS2doup(hObject, eventdata, handles)
AS2doup = get(hObject, 'Value');

handles.guidata.AS2doup = AS2doup;
guidata(hObject,handles)
%---
function AS1updt(hObject, eventdata, handles)
AS1updt = get(hObject, 'Value');

handles.guidata.AS1updt = AS1updt;
guidata(hObject,handles)
%---
function AS2updt(hObject, eventdata, handles)
AS2updt = get(hObject, 'Value');

handles.guidata.AS2updt = AS2updt;
guidata(hObject,handles)
%---
function FileATdataS1(hObject, eventdata, handles)
FileATdataS1 = get(hObject, 'Value');

handles.guidata.FileATdataS1 = FileATdataS1;
guidata(hObject,handles)
%---
function FileATdataS2(hObject, eventdata, handles)
FileATdataS2 = get(hObject, 'Value');

handles.guidata.FileATdataS2 = FileATdataS2;
guidata(hObject,handles)
%---
function ActS1scell(hObject, eventdata, handles)
ActS1scell = get(hObject, 'Value');

handles.guidata.ActS1scell = ActS1scell;
guidata(hObject,handles)
%---
function ActS2scell(hObject, eventdata, handles)
ActS2scell = get(hObject, 'Value');

handles.guidata.ActS2scell = ActS2scell;
guidata(hObject,handles)
%---
%=========================================================%
function MplotUp(hObject, eventdata, handles)
MplotUp = get(hObject, 'Value');

handles.guidata.MplotUp = MplotUp;
guidata(hObject,handles)
%---
function SplotUp(hObject, eventdata, handles)
SplotUp = get(hObject, 'Value');

handles.guidata.SplotUp = SplotUp;
guidata(hObject,handles)
%---
%=========================================================%
function KBrd(hObject, eventdata, handles)
KBrd = get(hObject, 'Value');

handles.guidata.KBrd = KBrd;
guidata(hObject,handles)
%---

