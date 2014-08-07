function varargout = BradsModelGUI(varargin)
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
%------------------------------------------------------------------------%
%========================================================================%


%  RunSimulation.
function RunSimulation_Callback(hObject, eventdata, handles)
% hObject    handle to RunSimulation (see GCBO)


dot1 = str2double(get(handles.GluR1dots,'String'));
dot2 = str2double(get(handles.GluR2dots,'String'));
dot3 = str2double(get(handles.Steps,'String'));
dot4 = str2double(get(handles.TimeStep,'String'));
dot5 = str2double(get(handles.Scale,'String'));
dot6 = str2double(get(handles.doLoops,'String'));

dot = [dot1 dot2 dot3 dot4 dot5 dot6];


dr1 = str2double(get(handles.esDGR1,'String'));
dr2 = str2double(get(handles.spineDGR1,'String'));
dr3 = str2double(get(handles.PERIDGR1,'String'));
dr4 = str2double(get(handles.PSDDGR1,'String'));
dr5 = str2double(get(handles.esDGR2,'String'));
dr6 = str2double(get(handles.spineDGR2,'String'));
dr7 = str2double(get(handles.PERIDGR2,'String'));
dr8 = str2double(get(handles.PSDDGR2,'String'));

dr = [dr1 dr2 dr3 dr4 dr5 dr6 dr7 dr8];


um1 = str2double(get(handles.denWidthX,'String'));
um2 = str2double(get(handles.denHeightY,'String'));
um3 = str2double(get(handles.PSD1um,'String'));
um4 = str2double(get(handles.PSD2um,'String'));
um5 = str2double(get(handles.PERI1um,'String'));
um6 = str2double(get(handles.PERI2um,'String'));

um = [um1 um2 um3 um4 um5 um6];


sap1 = str2double(get(handles.SAPdotsPSD1,'String'));
sap2 = str2double(get(handles.SAPdotsPSD2,'String'));
% sap3 = str2double(get(handles.SAPbetaPSD1,'String'));
% sap4 = str2double(get(handles.SAPtauPSD1,'String'));
% sap5 = str2double(get(handles.SAPL1PSD1,'String'));
% sap6 = str2double(get(handles.SAPbetaPSD2,'String'));
% sap7 = str2double(get(handles.SAPtauPSD2,'String'));
% sap8 = str2double(get(handles.SAPL1PSD2,'String'));
% sap9 = str2double(get(handles.SAPmuPSD1,'String'));
% sap10 = str2double(get(handles.SAPmuPSD2,'String'));
sap11 = str2double(get(handles.SAPdTPSD1,'String'));
sap12 = str2double(get(handles.SAPdTPSD2,'String'));
% sap13 = str2double(get(handles.SAPrhoPSD1,'String'));
% sap14 = str2double(get(handles.SAPrhoPSD2,'String'));
% sap15 = str2double(get(handles.SAPrPSD1,'String'));
% sap16 = str2double(get(handles.SAPrPSD2,'String'));
sap17 = get(handles.doDynamicLeP1,'Value');
sap18 = get(handles.doDynamicLeP2,'Value');
sap19 = get(handles.doLTPS1,'Value');
sap20 = get(handles.doLTPS2,'Value');

sap21 = str2double(get(handles.LonS1,'String'));
sap22 = str2double(get(handles.QonS1,'String'));
sap23 = str2double(get(handles.RonS1,'String'));
sap24 = str2double(get(handles.LoffS1,'String'));
sap25 = str2double(get(handles.QoffS1,'String'));
sap26 = str2double(get(handles.RoffS1,'String'));
sap27 = str2double(get(handles.LonS2,'String'));
sap28 = str2double(get(handles.QonS2,'String'));
sap29 = str2double(get(handles.RonS2,'String'));
sap30 = str2double(get(handles.LoffS2,'String'));
sap31 = str2double(get(handles.QoffS2,'String'));
sap32 = str2double(get(handles.RoffS2,'String'));

% sap = [sap1 sap2 sap3 sap4 sap5 sap6 sap7 sap8 sap9 sap10...
% 	   sap11 sap12 sap13 sap14 sap15 sap16 sap17 sap18 sap19 sap20...
% 	   sap21 sap22 sap23 sap24 sap25 sap26 sap27 sap28 sap29 sap30 sap31 sap32];

sap = [sap1 sap2 0 0 0 0 0 0 0 0 ...
	   sap11 sap12 0 0 0 0 sap17 sap18 sap19 sap20...
	   sap21 sap22 sap23 sap24 sap25 sap26 sap27 sap28 sap29 sap30 sap31 sap32];

   
   
hr1 = str2double(get(handles.HomeostaticLo,'String'));
hr2 = str2double(get(handles.HomeostaticHi,'String'));

hr = [hr1 hr2];


doUse1 = get(handles.useGluR1,'Value');
doUse2 = get(handles.useGluR2,'Value');
doUse3 = get(handles.runSAPPSD1,'Value');
doUse4 = get(handles.runSAPPSD2,'Value');

doUse = [doUse1 doUse2 doUse3 doUse4];


doRun1 = get(handles.run2Dplot,'Value');
doRun2 = get(handles.run3Dplot,'Value');
doRun3 = get(handles.runMSDtest,'Value');
doRun4 = get(handles.runTraceSingleDot,'Value');
doRun5 = get(handles.runMSDpopup,'Value');		% POPUP
doRun6 = get(handles.runHomeostatic,'Value');
doRun7 = get(handles.runPoissonsBox,'Value');
doRun8 = get(handles.doFieldFig,'Value');
doRun9 = get(handles.doSlotColormap,'Value');
doRun10 = get(handles.doSpike,'Value');
doRun11 = get(handles.doProfile,'Value');

doRun = [doRun1 doRun2 doRun3 doRun4 doRun5 doRun6 doRun7 ...
	     doRun8 doRun9 doRun10 doRun11];

box1 = str2double(get(handles.GraphTime,'String'));
box2 = str2double(get(handles.AllowedTime,'String'));
box3 = str2double(get(handles.NumLoops,'String'));
box4 = get(handles.SteadyState,'Value');

box = [box1 box2 box3 box4];



%{
slt1 = str2double(get(handles.G1PSDslotN,'String'));
slt2 = str2double(get(handles.G1PERIslotN,'String'));
slt3 = str2double(get(handles.G2PSDslotN,'String'));
slt4 = str2double(get(handles.G2PERIslotN,'String'));

slt = [slt1 slt2 slt3 slt4];



ko1 = str2double(get(handles.KonSpi1PSDGR2,'String'));
ko2 = str2double(get(handles.KoffSpi1PSDGR2,'String'));
ko3 = str2double(get(handles.KonSpi1PERIGR2,'String'));
ko4 = str2double(get(handles.KoffSpi1PERIGR2,'String'));
ko5 = str2double(get(handles.KonSpi2PSDGR2,'String'));
ko6 = str2double(get(handles.KoffSpi2PSDGR2,'String'));
ko7 = str2double(get(handles.KonSpi2PERIGR2,'String'));
ko8 = str2double(get(handles.KoffSpi2PERIGR2,'String'));

ko9 = str2double(get(handles.KonSpi1PSDGR1,'String'));
ko10 = str2double(get(handles.KoffSpi1PSDGR1,'String'));
ko11 = str2double(get(handles.KonSpi1PERIGR1,'String'));
ko12 = str2double(get(handles.KoffSpi1PERIGR1,'String'));
ko13 = str2double(get(handles.KonSpi2PSDGR1,'String'));
ko14 = str2double(get(handles.KoffSpi2PSDGR1,'String'));
ko15 = str2double(get(handles.KonSpi2PERIGR1,'String'));
ko16 = str2double(get(handles.KoffSpi2PERIGR1,'String'));

ko = [ko1 ko2 ko3 ko4 ko5 ko6 ko7 ko8 ko9 ko10 ko11 ko12 ko13 ko14 ko15 ko16];
%}
ko = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
slt = [1 1 1 1];



doKo1 = get(handles.useGluR1slots,'Value');
doKo2 = get(handles.useGluR2slots,'Value');
doKo3 = get(handles.useSlotsPSDGR1,'Value');
doKo4 = get(handles.useSlotsPSDGR2,'Value');
doKo5 = get(handles.useSlotsPERIGR1,'Value');
doKo6 = get(handles.useSlotsPERIGR2,'Value');

doKo = [doKo1 doKo2 doKo3 doKo4 doKo5 doKo6];


GT1 = str2double(get(handles.GT1on,'String'));
GT2 = str2double(get(handles.GT1off,'String'));
GT3 = str2double(get(handles.GT1LTPv,'String'));
GT4 = str2double(get(handles.GT2on,'String'));
GT5 = str2double(get(handles.GT2off,'String'));
GT6 = str2double(get(handles.GT2LTPv,'String'));

G1Tab = get(handles.GT1masktab,'Data');
G2Tab = get(handles.GT2masktab,'Data');

GT = [GT1 GT2 GT3 GT4 GT5 GT6];
GTab = {G1Tab, G2Tab};


stky1 = str2double(get(handles.G1STBASE,'String'));
stky2 = 10; %str2double(get(handles.G1RTBASE,'String'));
stky3 = str2double(get(handles.G1STLTP,'String'));
stky4 = 10; %str2double(get(handles.G1RTLTP,'String'));
stky5 = str2double(get(handles.G1BSMu,'String'));
stky6 = str2double(get(handles.G1LSMu,'String'));

stky7 = str2double(get(handles.G2STBASE,'String'));
stky8 = 10; %str2double(get(handles.G2RTBASE,'String'));
stky9 = str2double(get(handles.G2STLTP,'String'));
stky10 = 10;% str2double(get(handles.G2RTLTP,'String'));
stky11 = str2double(get(handles.G2BSMu,'String'));
stky12 = str2double(get(handles.G2LSMu,'String'));

stky13 = str2double(get(handles.G1TW1,'String'));
stky14 = str2double(get(handles.G1TW2,'String'));
stky15 = str2double(get(handles.G1TW3,'String'));
stky16 = str2double(get(handles.G1TW4,'String'));
stky17 = str2double(get(handles.G2TW1,'String'));
stky18 = str2double(get(handles.G2TW2,'String'));
stky19 = str2double(get(handles.G2TW3,'String'));
stky20 = str2double(get(handles.G2TW4,'String'));
stky21 = str2double(get(handles.SAPPADPSD1,'String'));
stky22 = str2double(get(handles.SAPPADPSD2,'String'));


stky = [stky1 stky2 stky3 stky4 stky5 stky6 stky7 stky8 stky9...
	    stky10 stky11 stky12 stky13 stky14 stky15 stky16 stky17...
		stky18 stky19 stky20 stky21 stky22];


% handles.output = ReDiClus(dot,dr,um,sap,hr,ko,doUse,doRun,doKo,box,slt,stky);
handles.output = MAINBOX(dot,dr,um,sap,hr,ko,doUse,doRun,doKo,box,slt,stky,GT,GTab);


%  RunPoisBox.
function RunPoisBox_Callback(hObject, eventdata, handles)
% hObject    handle to RunPoisBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% mass = handles.metricdata.density * handles.metricdata.volume;



num1 = str2double(get(handles.GluR2dots,'String'));
num2 = str2double(get(handles.GluR1dots,'String'));
num3 = str2double(get(handles.GraphTime,'String'));
num4 = str2double(get(handles.AllowedTime,'String'));
num5 = str2double(get(handles.spiD,'String'));
num6 = str2double(get(handles.psd2D,'String'));
num7 = str2double(get(handles.psd1D,'String'));
num8 = str2double(get(handles.NumLoops,'String'));
nums = [num1 num2 num3 num4 num5 num6 num7 num8];

dothis1 = get(handles.SteadyState,'Value');
dothis2 = get(handles.doLiveSim,'Value');
dothis = [dothis1 dothis2];

handles.output = PoissonsBox(nums,dothis);


%========================================================================%
%------------------------------------------------------------------------%
%========================================================================%
%------------------------------------------------------------------------%
%========================================================================%


% --- Executes during object creation, after setting all properties.
function GluR2dots_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GluR2dots (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function GluR2dots_Callback(hObject, eventdata, handles)

GluR2dots = str2double(get(hObject, 'String'));
if isnan(GluR2dots)
    set(hObject, 'String', 80);
    errordlg('Input must be a number','Error');
end

% Save the new GluR2dots value
handles.guidata.GluR2dots = GluR2dots;
guidata(hObject,handles)



function GraphTime_Callback(hObject, eventdata, handles)
GraphTime = str2double(get(hObject, 'String'));
if isnan(GraphTime)
    set(hObject, 'String', 30);
    errordlg('Input must be a number','Error');
end

% Save the new GraphTime value
handles.guidata.GraphTime = GraphTime;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function GraphTime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GraphTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function AllowedTime_Callback(hObject, eventdata, handles)
AllowedTime = str2double(get(hObject, 'String'));
if isnan(AllowedTime)
    set(hObject, 'String', 200);
    errordlg('Input must be a number','Error');
end

% Save the new AllowedTime value
handles.guidata.AllowedTime = AllowedTime;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function AllowedTime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AllowedTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function GluR1dots_Callback(hObject, eventdata, handles)
GluR1dots = str2double(get(hObject, 'String'));
if isnan(GluR1dots)
    set(hObject, 'String', 20);
    errordlg('Input must be a number','Error');
end

% Save the new GluR1dots value
handles.guidata.GluR1dots = GluR1dots;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function GluR1dots_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GluR1dots (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function spiD_Callback(hObject, eventdata, handles)
spiD = str2double(get(hObject, 'String'));
if isnan(spiD)
    set(hObject, 'String', 0.1);
    errordlg('Input must be a number','Error');
end

% Save the new spiD value
handles.guidata.spiD = spiD;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function spiD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to spiD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function psd2D_Callback(hObject, eventdata, handles)
psd2D = str2double(get(hObject, 'String'));
if isnan(psd2D)
    set(hObject, 'String', 0.01);
    errordlg('Input must be a number','Error');
end

% Save the new psd2D value
handles.guidata.psd2D = psd2D;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function psd2D_CreateFcn(hObject, eventdata, handles)
% hObject    handle to psd2D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function psd1D_Callback(hObject, eventdata, handles)
psd1D = str2double(get(hObject, 'String'));
if isnan(psd1D)
    set(hObject, 'String', 0.001);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.psd1D = psd1D;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function psd1D_CreateFcn(hObject, eventdata, handles)
% hObject    handle to psd1D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function NumLoops_Callback(hObject, eventdata, handles)
NumLoops = str2double(get(hObject, 'String'));
if isnan(NumLoops)
    set(hObject, 'String', 5);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.NumLoops = NumLoops;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function NumLoops_CreateFcn(hObject, eventdata, handles)
% hObject    handle to psd1D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%  SteadyState.
function SteadyState_Callback(hObject, eventdata, handles)
% Hint: get(hObject,'Value') returns toggle state of SteadyState

SteadyState = get(hObject, 'Value');

% Save the new SteadyState value
handles.guidata.SteadyState = SteadyState;
guidata(hObject,handles)


%  doLiveSim.
function doLiveSim_Callback(hObject, eventdata, handles)
% Hint: get(hObject,'Value') returns toggle state of doLiveSim

doLiveSim = get(hObject, 'Value');

% Save the new doLiveSim value
handles.guidata.doLiveSim = doLiveSim;
guidata(hObject,handles)


%=========================================================%

function doLoops_Callback(hObject, eventdata, handles)
doLoops = str2double(get(hObject, 'String'));
if isnan(doLoops)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.doLoops = doLoops;
guidata(hObject,handles)


function Steps_Callback(hObject, eventdata, handles)
Steps = str2double(get(hObject, 'String'));
if isnan(Steps)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.Steps = Steps;
guidata(hObject,handles)


function TimeStep_Callback(hObject, eventdata, handles)
TimeStep = str2double(get(hObject, 'String'));
if isnan(TimeStep)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.TimeStep = TimeStep;
guidata(hObject,handles)


function Scale_Callback(hObject, eventdata, handles)
Scale = str2double(get(hObject, 'String'));
if isnan(Scale)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.Scale = Scale;
guidata(hObject,handles)


%=========================================================%

%  runMSDtest.
function runMSDtest_Callback(hObject, eventdata, handles)
runMSDtest = get(hObject, 'Value');

% Save the new SteadyState value
handles.guidata.runMSDtest = runMSDtest;
guidata(hObject,handles)



%  runTraceSingleDot.
function runTraceSingleDot_Callback(hObject, eventdata, handles)
runTraceSingleDot = get(hObject, 'Value');

% Save the new SteadyState value
handles.guidata.runTraceSingleDot = runTraceSingleDot;
guidata(hObject,handles)


% --- Executes on selection change in runMSDpopup.
function runMSDpopup_Callback(hObject, eventdata, handles)
runMSDpopup = get(hObject, 'Value');

% Save the new SteadyState value
handles.guidata.runMSDpopup = runMSDpopup;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function runMSDpopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to runMSDpopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%  runHomeostatic.
function runHomeostatic_Callback(hObject, eventdata, handles)
runHomeostatic = get(hObject, 'Value');

% Save the new SteadyState value
handles.guidata.runHomeostatic = runHomeostatic;
guidata(hObject,handles)



function HomeostaticHi_Callback(hObject, eventdata, handles)
HomeostaticHi = str2double(get(hObject, 'String'));
if isnan(HomeostaticHi)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.HomeostaticHi = HomeostaticHi;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function HomeostaticHi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to HomeostaticHi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function HomeostaticLo_Callback(hObject, eventdata, handles)
HomeostaticLo = str2double(get(hObject, 'String'));
if isnan(HomeostaticLo)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.HomeostaticLo = HomeostaticLo;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function HomeostaticLo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to HomeostaticLo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%  run2Dplot.
function run2Dplot_Callback(hObject, eventdata, handles)
run2Dplot = get(hObject, 'Value');

% Save the new SteadyState value
handles.guidata.run2Dplot = run2Dplot;
guidata(hObject,handles)


%  run3Dplot.
function run3Dplot_Callback(hObject, eventdata, handles)
run3Dplot = get(hObject, 'Value');

% Save the new SteadyState value
handles.guidata.run3Dplot = run3Dplot;
guidata(hObject,handles)





% --- Executes during object creation, after setting all properties.
function SAPbetaPSD1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SAPbetaPSD1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





% --- Executes during object creation, after setting all properties.
function SAPtauPSD1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SAPtauPSD1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





% --- Executes during object creation, after setting all properties.
function SAPL1PSD1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SAPL1PSD1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function SAPdotsPSD1_Callback(hObject, eventdata, handles)
SAPdotsPSD1 = str2double(get(hObject, 'String'));
if isnan(SAPdotsPSD1)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.SAPdotsPSD1 = SAPdotsPSD1;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function SAPdotsPSD1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SAPdotsPSD1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





% --- Executes during object creation, after setting all properties.
function SAPbetaPSD2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SAPbetaPSD2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





% --- Executes during object creation, after setting all properties.
function SAPtauPSD2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SAPtauPSD2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





% --- Executes during object creation, after setting all properties.
function SAPL1PSD2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SAPL1PSD2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function SAPdotsPSD2_Callback(hObject, eventdata, handles)
SAPdotsPSD2 = str2double(get(hObject, 'String'));
if isnan(SAPdotsPSD2)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.SAPdotsPSD2 = SAPdotsPSD2;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function SAPdotsPSD2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SAPdotsPSD2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%  runSAPPSD1.
function runSAPPSD1_Callback(hObject, eventdata, handles)
runSAPPSD1 = get(hObject, 'Value');

% Save the new SteadyState value
handles.guidata.runSAPPSD1 = runSAPPSD1;
guidata(hObject,handles)


%  runSAPPSD2.
function runSAPPSD2_Callback(hObject, eventdata, handles)
runSAPPSD2 = get(hObject, 'Value');

% Save the new SteadyState value
handles.guidata.runSAPPSD2 = runSAPPSD2;
guidata(hObject,handles)


%  runSAPPSD1.
function doDynamicLeP1_Callback(hObject, eventdata, handles)
doDynamicLeP1 = get(hObject, 'Value');

% Save the new SteadyState value
handles.guidata.doDynamicLeP1 = doDynamicLeP1;
guidata(hObject,handles)


%  runSAPPSD2.
function doDynamicLeP2_Callback(hObject, eventdata, handles)
doDynamicLeP2 = get(hObject, 'Value');

% Save the new SteadyState value
handles.guidata.doDynamicLeP2 = doDynamicLeP2;
guidata(hObject,handles)


%  runSAPPSD1.
function doLTPS1_Callback(hObject, eventdata, handles)
doLTPS1 = get(hObject, 'Value');

% Save the new SteadyState value
handles.guidata.doLTPS1 = doLTPS1;
guidata(hObject,handles)


%  runSAPPSD2.
function doLTPS2_Callback(hObject, eventdata, handles)
doLTPS2 = get(hObject, 'Value');

% Save the new SteadyState value
handles.guidata.doLTPS2 = doLTPS2;
guidata(hObject,handles)


function denWidthX_Callback(hObject, eventdata, handles)
denWidthX = str2double(get(hObject, 'String'));
if isnan(denWidthX)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.denWidthX = denWidthX;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function denWidthX_CreateFcn(hObject, eventdata, handles)
% hObject    handle to denWidthX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function denHeightY_Callback(hObject, eventdata, handles)
denHeightY = str2double(get(hObject, 'String'));
if isnan(denHeightY)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.denHeightY = denHeightY;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function denHeightY_CreateFcn(hObject, eventdata, handles)
% hObject    handle to denHeightY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function PSD1um_Callback(hObject, eventdata, handles)
PSD1um = str2double(get(hObject, 'String'));
if isnan(PSD1um)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.PSD1um = PSD1um;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function PSD1um_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PSD1um (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function PSD2um_Callback(hObject, eventdata, handles)
PSD2um = str2double(get(hObject, 'String'));
if isnan(PSD2um)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.PSD2um = PSD2um;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function PSD2um_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PSD2um (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function PERI1um_Callback(hObject, eventdata, handles)
PERI1um = str2double(get(hObject, 'String'));
if isnan(PERI1um)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.PERI1um = PERI1um;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function PERI1um_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PERI1um (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function PERI2um_Callback(hObject, eventdata, handles)
PERI2um = str2double(get(hObject, 'String'));
if isnan(PERI2um)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.PERI2um = PERI2um;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function PERI2um_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PERI2um (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function spineDGR1_Callback(hObject, eventdata, handles)
spineDGR1 = str2double(get(hObject, 'String'));
if isnan(spineDGR1)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.spineDGR1 = spineDGR1;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function spineDGR1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to spineDGR1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function PERIDGR1_Callback(hObject, eventdata, handles)
PERIDGR1 = str2double(get(hObject, 'String'));
if isnan(PERIDGR1)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.PERIDGR1 = PERIDGR1;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function PERIDGR1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PERIDGR1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function PSDDGR1_Callback(hObject, eventdata, handles)
PSDDGR1 = str2double(get(hObject, 'String'));
if isnan(PSDDGR1)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.PSDDGR1 = PSDDGR1;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function PSDDGR1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PSDDGR1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function esDGR1_Callback(hObject, eventdata, handles)
esDGR1 = str2double(get(hObject, 'String'));
if isnan(esDGR1)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.esDGR1 = esDGR1;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function esDGR1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to esDGR1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function spineDGR2_Callback(hObject, eventdata, handles)
spineDGR2 = str2double(get(hObject, 'String'));
if isnan(spineDGR2)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.spineDGR2 = spineDGR2;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function spineDGR2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to spineDGR2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function PERIDGR2_Callback(hObject, eventdata, handles)
PERIDGR2 = str2double(get(hObject, 'String'));
if isnan(PERIDGR2)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.PERIDGR2 = PERIDGR2;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function PERIDGR2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PERIDGR2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function PSDDGR2_Callback(hObject, eventdata, handles)
PSDDGR2 = str2double(get(hObject, 'String'));
if isnan(PSDDGR2)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.PSDDGR2 = PSDDGR2;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function PSDDGR2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PSDDGR2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function esDGR2_Callback(hObject, eventdata, handles)
esDGR2 = str2double(get(hObject, 'String'));
if isnan(esDGR2)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.esDGR2 = esDGR2;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function esDGR2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to esDGR2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%  runPoissonsBox.
function runPoissonsBox_Callback(hObject, eventdata, handles)
runPoissonsBox = get(hObject, 'Value');

% Save the new SteadyState value
handles.guidata.runPoissonsBox = runPoissonsBox;
guidata(hObject,handles)


%  useGluR1.
function useGluR1_Callback(hObject, eventdata, handles)
useGluR1 = get(hObject, 'Value');

% Save the new SteadyState value
handles.guidata.useGluR1 = useGluR1;
guidata(hObject,handles)


%  useGluR2.
function useGluR2_Callback(hObject, eventdata, handles)
useGluR2 = get(hObject, 'Value');

% Save the new SteadyState value
handles.guidata.useGluR2 = useGluR2;
guidata(hObject,handles)



function KonSpi1PSDGR2_Callback(hObject, eventdata, handles)
KonSpi1PSDGR2 = str2double(get(hObject, 'String'));
if isnan(KonSpi1PSDGR2)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.KonSpi1PSDGR2 = KonSpi1PSDGR2;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function KonSpi1PSDGR2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to KonSpi1PSDGR2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function KoffSpi1PSDGR2_Callback(hObject, eventdata, handles)
KoffSpi1PSDGR2 = str2double(get(hObject, 'String'));
if isnan(KoffSpi1PSDGR2)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.KoffSpi1PSDGR2 = KoffSpi1PSDGR2;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function KoffSpi1PSDGR2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to KoffSpi1PSDGR2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function KonSpi1PERIGR2_Callback(hObject, eventdata, handles)
KonSpi1PERIGR2 = str2double(get(hObject, 'String'));
if isnan(KonSpi1PERIGR2)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.KonSpi1PERIGR2 = KonSpi1PERIGR2;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function KonSpi1PERIGR2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to KonSpi1PERIGR2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function KoffSpi1PERIGR2_Callback(hObject, eventdata, handles)
KoffSpi1PERIGR2 = str2double(get(hObject, 'String'));
if isnan(KoffSpi1PERIGR2)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.KoffSpi1PERIGR2 = KoffSpi1PERIGR2;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function KoffSpi1PERIGR2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to KoffSpi1PERIGR2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function KonSpi2PSDGR2_Callback(hObject, eventdata, handles)
KonSpi2PSDGR2 = str2double(get(hObject, 'String'));
if isnan(KonSpi2PSDGR2)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.KonSpi2PSDGR2 = KonSpi2PSDGR2;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function KonSpi2PSDGR2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to KonSpi2PSDGR2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function KoffSpi2PSDGR2_Callback(hObject, eventdata, handles)
KoffSpi2PSDGR2 = str2double(get(hObject, 'String'));
if isnan(KoffSpi2PSDGR2)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.KoffSpi2PSDGR2 = KoffSpi2PSDGR2;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function KoffSpi2PSDGR2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to KoffSpi2PSDGR2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function KonSpi2PERIGR2_Callback(hObject, eventdata, handles)
KonSpi2PERIGR2 = str2double(get(hObject, 'String'));
if isnan(KonSpi2PERIGR2)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.KonSpi2PERIGR2 = KonSpi2PERIGR2;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function KonSpi2PERIGR2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to KonSpi2PERIGR2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function KoffSpi2PERIGR2_Callback(hObject, eventdata, handles)
KoffSpi2PERIGR2 = str2double(get(hObject, 'String'));
if isnan(KoffSpi2PERIGR2)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.KoffSpi2PERIGR2 = KoffSpi2PERIGR2;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function KoffSpi2PERIGR2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to KoffSpi2PERIGR2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function KonSpi1PSDGR1_Callback(hObject, eventdata, handles)
KonSpi1PSDGR1 = str2double(get(hObject, 'String'));
if isnan(KonSpi1PSDGR1)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.KonSpi1PSDGR1 = KonSpi1PSDGR1;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function KonSpi1PSDGR1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to KonSpi1PSDGR1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function KoffSpi1PSDGR1_Callback(hObject, eventdata, handles)
KoffSpi1PSDGR1 = str2double(get(hObject, 'String'));
if isnan(KoffSpi1PSDGR1)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.KoffSpi1PSDGR1 = KoffSpi1PSDGR1;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function KoffSpi1PSDGR1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to KoffSpi1PSDGR1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function KonSpi1PERIGR1_Callback(hObject, eventdata, handles)
KonSpi1PERIGR1 = str2double(get(hObject, 'String'));
if isnan(KonSpi1PERIGR1)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.KonSpi1PERIGR1 = KonSpi1PERIGR1;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function KonSpi1PERIGR1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to KonSpi1PERIGR1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function KoffSpi1PERIGR1_Callback(hObject, eventdata, handles)
KoffSpi1PERIGR1 = str2double(get(hObject, 'String'));
if isnan(KoffSpi1PERIGR1)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.KoffSpi1PERIGR1 = KoffSpi1PERIGR1;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function KoffSpi1PERIGR1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to KoffSpi1PERIGR1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function KonSpi2PSDGR1_Callback(hObject, eventdata, handles)
KonSpi2PSDGR1 = str2double(get(hObject, 'String'));
if isnan(KonSpi2PSDGR1)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.KonSpi2PSDGR1 = KonSpi2PSDGR1;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function KonSpi2PSDGR1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to KonSpi2PSDGR1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function KoffSpi2PSDGR1_Callback(hObject, eventdata, handles)
KoffSpi2PSDGR1 = str2double(get(hObject, 'String'));
if isnan(KoffSpi2PSDGR1)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.KoffSpi2PSDGR1 = KoffSpi2PSDGR1;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function KoffSpi2PSDGR1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to KoffSpi2PSDGR1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function KonSpi2PERIGR1_Callback(hObject, eventdata, handles)
KonSpi2PERIGR1 = str2double(get(hObject, 'String'));
if isnan(KonSpi2PERIGR1)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.KonSpi2PERIGR1 = KonSpi2PERIGR1;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function KonSpi2PERIGR1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to KonSpi2PERIGR1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function KoffSpi2PERIGR1_Callback(hObject, eventdata, handles)
KoffSpi2PERIGR1 = str2double(get(hObject, 'String'));
if isnan(KoffSpi2PERIGR1)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.KoffSpi2PERIGR1 = KoffSpi2PERIGR1;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function KoffSpi2PERIGR1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to KoffSpi2PERIGR1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%  useGluR1slots.
function useGluR1slots_Callback(hObject, eventdata, handles)
useGluR1slots = get(hObject, 'Value');

% Save the new SteadyState value
handles.guidata.useGluR1slots = useGluR1slots;
guidata(hObject,handles)


%  useSlotsPSD1GR1.
function useSlotsPSDGR1_Callback(hObject, eventdata, handles)
useSlotsPSDGR1 = get(hObject, 'Value');

% Save the new SteadyState value
handles.guidata.useSlotsPSDGR1 = useSlotsPSDGR1;
guidata(hObject,handles)


%  useSlotsPERI1GR1.
function useSlotsPERIGR1_Callback(hObject, eventdata, handles)
useSlotsPERIGR1 = get(hObject, 'Value');

% Save the new SteadyState value
handles.guidata.useSlotsPERIGR1 = useSlotsPERIGR1;
guidata(hObject,handles)


%  useGluR2slots.
function useGluR2slots_Callback(hObject, eventdata, handles)
useGluR2slots = get(hObject, 'Value');

% Save the new SteadyState value
handles.guidata.useGluR2slots = useGluR2slots;
guidata(hObject,handles)


%  useSlotsPSD1GR2.
function useSlotsPSDGR2_Callback(hObject, eventdata, handles)
useSlotsPSDGR2 = get(hObject, 'Value');

% Save the new SteadyState value
handles.guidata.useSlotsPSDGR2 = useSlotsPSDGR2;
guidata(hObject,handles)


%  useSlotsPERI1GR2.
function useSlotsPERIGR2_Callback(hObject, eventdata, handles)
useSlotsPERIGR2 = get(hObject, 'Value');

% Save the new SteadyState value
handles.guidata.useSlotsPERIGR2 = useSlotsPERIGR2;
guidata(hObject,handles)

%  useSlotsPERI1GR2.
function SAPrunSLOTS_Callback(hObject, eventdata, handles)
SAPrunSLOTS = get(hObject, 'Value');

% Save the new SteadyState value
handles.guidata.SAPrunSLOTS = SAPrunSLOTS;
guidata(hObject,handles)

%=========================================================%

function G1PSDslotN_Callback(hObject, eventdata, handles)
G1PSDslotN = str2double(get(hObject, 'String'));
if isnan(G1PSDslotN)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.G1PSDslotN = G1PSDslotN;
guidata(hObject,handles)

function G1PERIslotN_Callback(hObject, eventdata, handles)
G1PERIslotN = str2double(get(hObject, 'String'));
if isnan(G1PERIslotN)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.G1PERIslotN = G1PERIslotN;
guidata(hObject,handles)

function G2PSDslotN_Callback(hObject, eventdata, handles)
G2PSDslotN = str2double(get(hObject, 'String'));
if isnan(G2PSDslotN)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.G2PSDslotN = G2PSDslotN;
guidata(hObject,handles)

function G2PERIslotN_Callback(hObject, eventdata, handles)
G2PERIslotN = str2double(get(hObject, 'String'));
if isnan(G2PERIslotN)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.G2PERIslotN = G2PERIslotN;
guidata(hObject,handles)

%=========================================================%

function G1TW1_Callback(hObject, eventdata, handles)
G1TW1 = str2double(get(hObject, 'String'));
if isnan(G1TW1)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.G1TW1 = G1TW1;
guidata(hObject,handles)

function G1TW2_Callback(hObject, eventdata, handles)
G1TW2 = str2double(get(hObject, 'String'));
if isnan(G1TW2)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.G1TW2 = G1TW2;
guidata(hObject,handles)

function G1TW3_Callback(hObject, eventdata, handles)
G1TW3 = str2double(get(hObject, 'String'));
if isnan(G1TW3)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.G1TW3 = G1TW3;
guidata(hObject,handles)

function G1TW4_Callback(hObject, eventdata, handles)
G1TW4 = str2double(get(hObject, 'String'));
if isnan(G1TW4)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.G1TW4 = G1TW4;
guidata(hObject,handles)

function G2TW1_Callback(hObject, eventdata, handles)
G2TW1 = str2double(get(hObject, 'String'));
if isnan(G2TW1)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.G2TW1 = G2TW1;
guidata(hObject,handles)

function G2TW2_Callback(hObject, eventdata, handles)
G2TW2 = str2double(get(hObject, 'String'));
if isnan(G2TW2)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.G2TW2 = G2TW2;
guidata(hObject,handles)

function G2TW3_Callback(hObject, eventdata, handles)
G2TW3 = str2double(get(hObject, 'String'));
if isnan(G2TW3)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.G2TW3 = G2TW3;
guidata(hObject,handles)

function G2TW4_Callback(hObject, eventdata, handles)
G2TW4 = str2double(get(hObject, 'String'));
if isnan(G2TW4)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.G2TW4 = G2TW4;
guidata(hObject,handles)

function SAPPADPSD1_Callback(hObject, eventdata, handles)
SAPPADPSD1 = str2double(get(hObject, 'String'));
if isnan(SAPPADPSD1)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.SAPPADPSD1 = SAPPADPSD1;
guidata(hObject,handles)

function SAPPADPSD2_Callback(hObject, eventdata, handles)
SAPPADPSD2 = str2double(get(hObject, 'String'));
if isnan(SAPPADPSD2)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.SAPPADPSD2 = SAPPADPSD2;
guidata(hObject,handles)

%=========================================================%

function G1STBASE_Callback(hObject, eventdata, handles)
G1STKY = str2double(get(hObject, 'String'));
if isnan(G1STKY)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.G1STKY = G1STKY;
guidata(hObject,handles)

function G1RTBASE_Callback(hObject, eventdata, handles)
G1STAB = str2double(get(hObject, 'String'));
if isnan(G1STAB)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.G1STAB = G1STAB;
guidata(hObject,handles)

function G1STLTP_Callback(hObject, eventdata, handles)
G1STLTP = str2double(get(hObject, 'String'));
if isnan(G1STLTP)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.G1STLTP = G1STLTP;
guidata(hObject,handles)

function G1RTLTP_Callback(hObject, eventdata, handles)
G1RTLTP = str2double(get(hObject, 'String'));
if isnan(G1RTLTP)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.G1RTLTP = G1RTLTP;
guidata(hObject,handles)

function G2STBASE_Callback(hObject, eventdata, handles)
G2STBASE = str2double(get(hObject, 'String'));
if isnan(G2STBASE)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.G2STBASE = G2STBASE;
guidata(hObject,handles)

function G2RTBASE_Callback(hObject, eventdata, handles)
G2RTBASE = str2double(get(hObject, 'String'));
if isnan(G2RTBASE)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.G2RTBASE = G2RTBASE;
guidata(hObject,handles)

function G2STLTP_Callback(hObject, eventdata, handles)
G2STLTP = str2double(get(hObject, 'String'));
if isnan(G2STLTP)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.G2STLTP = G2STLTP;
guidata(hObject,handles)

function G2RTLTP_Callback(hObject, eventdata, handles)
G2RTLTP = str2double(get(hObject, 'String'));
if isnan(G2RTLTP)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.G2RTLTP = G2RTLTP;
guidata(hObject,handles)

function G1BSMu_Callback(hObject, eventdata, handles)
G1BSMu = str2double(get(hObject, 'String'));
if isnan(G1BSMu)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.G1BSMu = G1BSMu;
guidata(hObject,handles)

function G1LSMu_Callback(hObject, eventdata, handles)
G1LSMu = str2double(get(hObject, 'String'));
if isnan(G1LSMu)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.G1LSMu = G1LSMu;
guidata(hObject,handles)

function G2BSMu_Callback(hObject, eventdata, handles)
G2BSMu = str2double(get(hObject, 'String'));
if isnan(G2BSMu)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.G2BSMu = G2BSMu;
guidata(hObject,handles)

function G2LSMu_Callback(hObject, eventdata, handles)
G2LSMu = str2double(get(hObject, 'String'));
if isnan(G2LSMu)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new psd1D value
handles.guidata.G2LSMu = G2LSMu;
guidata(hObject,handles)


%=========================================================%
% --- Checkbox for doFieldFig
function doFieldFig_Callback(hObject, eventdata, handles)
doFieldFig = get(hObject, 'Value');

handles.guidata.doFieldFig = doFieldFig;
guidata(hObject,handles)
 


% --- Checkbox for doSlotColormap
function doSlotColormap_Callback(hObject, eventdata, handles)
doSlotColormap = get(hObject, 'Value');

handles.guidata.doSlotColormap = doSlotColormap;
guidata(hObject,handles)
%=========================================================%
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
function QonS1(hObject, eventdata, handles)
QonS1 = str2double(get(hObject, 'String'));
if isnan(QonS1)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

handles.guidata.QonS1 = QonS1;
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
function QoffS1(hObject, eventdata, handles)
QoffS1 = str2double(get(hObject, 'String'));
if isnan(QoffS1)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

handles.guidata.QoffS1 = QoffS1;
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
function QonS2(hObject, eventdata, handles)
QonS2 = str2double(get(hObject, 'String'));
if isnan(QonS2)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end
 
handles.guidata.QonS2 = QonS2;
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
function QoffS2(hObject, eventdata, handles)
QoffS2 = str2double(get(hObject, 'String'));
if isnan(QoffS2)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end
 
handles.guidata.QoffS2 = QoffS2;
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

function SAPdTPSD1_Callback(hObject, eventdata, handles)
SAPdTPSD1 = str2double(get(hObject, 'String'));
if isnan(SAPdTPSD1)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end
 
% Save the new psd1D value
handles.guidata.SAPdTPSD1 = SAPdTPSD1;
guidata(hObject,handles)


function SAPdTPSD2_Callback(hObject, eventdata, handles)
SAPdTPSD2 = str2double(get(hObject, 'String'));
if isnan(SAPdTPSD2)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end
 
% Save the new psd1D value
handles.guidata.SAPdTPSD2 = SAPdTPSD2;
guidata(hObject,handles)

 %=========================================================%
%=========================================================%
%{
function SAPmuPSD1_Callback(hObject, eventdata, handles)
SAPmuPSD1 = str2double(get(hObject, 'String'));
if isnan(SAPmuPSD1)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end
 
% Save the new psd1D value
handles.guidata.SAPmuPSD1 = SAPmuPSD1;
guidata(hObject,handles)


function SAPmuPSD2_Callback(hObject, eventdata, handles)
SAPmuPSD2 = str2double(get(hObject, 'String'));
if isnan(SAPmuPSD2)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end
 
% Save the new psd1D value
handles.guidata.SAPmuPSD2 = SAPmuPSD2;
guidata(hObject,handles)


function SAPdTPSD1_Callback(hObject, eventdata, handles)
SAPdTPSD1 = str2double(get(hObject, 'String'));
if isnan(SAPdTPSD1)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end
 
% Save the new psd1D value
handles.guidata.SAPdTPSD1 = SAPdTPSD1;
guidata(hObject,handles)


function SAPdTPSD2_Callback(hObject, eventdata, handles)
SAPdTPSD2 = str2double(get(hObject, 'String'));
if isnan(SAPdTPSD2)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end
 
% Save the new psd1D value
handles.guidata.SAPdTPSD2 = SAPdTPSD2;
guidata(hObject,handles)


function SAPrhoPSD1_Callback(hObject, eventdata, handles)
SAPrhoPSD1 = str2double(get(hObject, 'String'));
if isnan(SAPrhoPSD1)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end
 
% Save the new psd1D value
handles.guidata.SAPrhoPSD1 = SAPrhoPSD1;
guidata(hObject,handles)


function SAPrhoPSD2_Callback(hObject, eventdata, handles)
SAPrhoPSD2 = str2double(get(hObject, 'String'));
if isnan(SAPrhoPSD2)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end
 
% Save the new psd1D value
handles.guidata.SAPrhoPSD2 = SAPrhoPSD2;
guidata(hObject,handles)


function SAPrPSD1_Callback(hObject, eventdata, handles)
SAPrPSD1 = str2double(get(hObject, 'String'));
if isnan(SAPrPSD1)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end
 
% Save the new psd1D value
handles.guidata.SAPrPSD1 = SAPrPSD1;
guidata(hObject,handles)


function SAPrPSD2_Callback(hObject, eventdata, handles)
SAPrPSD2 = str2double(get(hObject, 'String'));
if isnan(SAPrPSD2)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end
 
% Save the new psd1D value
handles.guidata.SAPrPSD2 = SAPrPSD2;
guidata(hObject,handles)
%}
%{.
function SAPmuPSD1_Callback(hObject, eventdata, handles)
SAPmuPSD1 = str2double(get(hObject, 'String'));
if isnan(SAPmuPSD1)
    SAPmuPSD1 = 0.5;
end
 
% Save the new psd1D value
handles.guidata.SAPmuPSD1 = SAPmuPSD1;
guidata(hObject,handles)


function SAPmuPSD2_Callback(hObject, eventdata, handles)
SAPmuPSD2 = str2double(get(hObject, 'String'));
if isnan(SAPmuPSD2)
    SAPmuPSD2 = 0.5;
end
 
% Save the new psd1D value
handles.guidata.SAPmuPSD2 = SAPmuPSD2;
guidata(hObject,handles)


function SAPrhoPSD1_Callback(hObject, eventdata, handles)
SAPrhoPSD1 = str2double(get(hObject, 'String'));
if isnan(SAPrhoPSD1)
    SAPrhoPSD1 = 0.9;
end
 
% Save the new psd1D value
handles.guidata.SAPrhoPSD1 = SAPrhoPSD1;
guidata(hObject,handles)


function SAPrhoPSD2_Callback(hObject, eventdata, handles)
SAPrhoPSD2 = str2double(get(hObject, 'String'));
if isnan(SAPrhoPSD2)
    SAPrhoPSD2 = 0.9;
end
 
% Save the new psd1D value
handles.guidata.SAPrhoPSD2 = SAPrhoPSD2;
guidata(hObject,handles)


function SAPrPSD1_Callback(hObject, eventdata, handles)
SAPrPSD1 = str2double(get(hObject, 'String'));
if isnan(SAPrPSD1)
    SAPrPSD1 = 10;
end
 
% Save the new psd1D value
handles.guidata.SAPrPSD1 = SAPrPSD1;
guidata(hObject,handles)


function SAPrPSD2_Callback(hObject, eventdata, handles)
SAPrPSD2 = str2double(get(hObject, 'String'));
if isnan(SAPrPSD2)
    SAPrPSD2 = 10;
end
 
% Save the new psd1D value
handles.guidata.SAPrPSD2 = SAPrPSD2;
guidata(hObject,handles)


function SAPbetaPSD1_Callback(hObject, eventdata, handles)
SAPbetaPSD1 = str2double(get(hObject, 'String'));
if isnan(SAPbetaPSD1)
    SAPbetaPSD1 = 60;
end
 
% Save the new psd1D value
handles.guidata.SAPbetaPSD1 = SAPbetaPSD1;
guidata(hObject,handles)

function SAPbetaPSD2_Callback(hObject, eventdata, handles)
SAPbetaPSD2 = str2double(get(hObject, 'String'));
if isnan(SAPbetaPSD2)
    SAPbetaPSD2 = 60;
end
 
% Save the new psd1D value
handles.guidata.SAPbetaPSD2 = SAPbetaPSD2;
guidata(hObject,handles)


function SAPtauPSD1_Callback(hObject, eventdata, handles)
SAPtauPSD1 = str2double(get(hObject, 'String'));
if isnan(SAPtauPSD1)
    SAPtauPSD1 = 0.8;
end
 
% Save the new psd1D value
handles.guidata.SAPtauPSD1 = SAPtauPSD1;
guidata(hObject,handles)

function SAPtauPSD2_Callback(hObject, eventdata, handles)
SAPtauPSD2 = str2double(get(hObject, 'String'));
if isnan(SAPtauPSD2)
    SAPtauPSD2 = 0.8;
end
 
% Save the new psd1D value
handles.guidata.SAPtauPSD2 = SAPtauPSD2;
guidata(hObject,handles)


function SAPL1PSD1_Callback(hObject, eventdata, handles)
SAPL1PSD1 = str2double(get(hObject, 'String'));
if isnan(SAPL1PSD1)
    SAPL1PSD1 = 2.0;
end
 
% Save the new psd1D value
handles.guidata.SAPL1PSD1 = SAPL1PSD1;
guidata(hObject,handles)

function SAPL1PSD2_Callback(hObject, eventdata, handles)
SAPL1PSD2 = str2double(get(hObject, 'String'));
if isnan(SAPL1PSD2)
    SAPL1PSD2 = 2.0;
end
 
% Save the new psd1D value
handles.guidata.SAPL1PSD2 = SAPL1PSD2;
guidata(hObject,handles)

%}
%=========================================================%

function doSpike(hObject, eventdata, handles)
doSpike = get(hObject, 'Value');

handles.guidata.doSpike = doSpike;
guidata(hObject,handles)

function doProfile(hObject, eventdata, handles)
doProfile = get(hObject, 'Value');

handles.guidata.doProfile = doProfile;
guidata(hObject,handles)

%=========================================================%

function GT1on(hObject, eventdata, handles)
GT1on = str2double(get(hObject, 'String'));
if isnan(GT1on)
    GT1on = 10;
end
 
% Save the new psd1D value
handles.guidata.GT1on = GT1on;
guidata(hObject,handles)

function GT1off(hObject, eventdata, handles)
GT1off = str2double(get(hObject, 'String'));
if isnan(GT1off)
    GT1off = -5;
end
 
% Save the new psd1D value
handles.guidata.GT1off = GT1off;
guidata(hObject,handles)

function GT1LTPv(hObject, eventdata, handles)
GT1LTPv = str2double(get(hObject, 'String'));
if isnan(GT1LTPv)
    GT1LTPv = .001;
end
 
% Save the new psd1D value
handles.guidata.GT1LTPv = GT1LTPv;
guidata(hObject,handles)

function GT2on(hObject, eventdata, handles)
GT2on = str2double(get(hObject, 'String'));
if isnan(GT2on)
    GT1on = 10;
end
 
% Save the new psd1D value
handles.guidata.GT2on = GT2on;
guidata(hObject,handles)

function GT2off(hObject, eventdata, handles)
GT2off = str2double(get(hObject, 'String'));
if isnan(GT2off)
    GT2off = -5;
end
 
% Save the new psd1D value
handles.guidata.GT2off = GT2off;
guidata(hObject,handles)

function GT2LTPv(hObject, eventdata, handles)
GT2LTPv = str2double(get(hObject, 'String'));
if isnan(GT2LTPv)
    GT2LTPv = .001;
end
 
% Save the new psd1D value
handles.guidata.GT2LTPv = GT2LTPv;
guidata(hObject,handles)

function GT1masktab(hObject, eventdata, handles)
GT1masktab = get(hObject, 'Data');
if isnan(GT1masktab)
    GT1masktab = logical([0 1 0; 1 1 1; 0 1 0]);
end
handles.guidata.GT1masktab = GT1masktab;
guidata(hObject,handles)

function GT2masktab(hObject, eventdata, handles)
GT2masktab = get(hObject, 'Data');
if isnan(GT2masktab)
    GT2masktab = logical([0 1 0; 1 1 1; 0 1 0]);
end
handles.guidata.GT2masktab = GT2masktab;
guidata(hObject,handles)
%=========================================================%






