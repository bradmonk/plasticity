function varargout = SAPORMGUI(varargin)
% SAPORMGUI MATLAB code for SAPORMGUI.fig
%      SAPORMGUI, by itself, creates a new SAPORMGUI or raises the existing
%      singleton*.
%
%      H = SAPORMGUI returns the handle to a new SAPORMGUI or the handle to
%      the existing singleton*.
%
%      SAPORMGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SAPORMGUI.M with the given input arguments.
%
%      SAPORMGUI('Property','Value',...) creates a new SAPORMGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SAPORMGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SAPORMGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SAPORMGUI

% Last Modified by GUIDE v2.5 13-Feb-2014 18:46:44

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SAPORMGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @SAPORMGUI_OutputFcn, ...
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

% --- Executes just before SAPORMGUI is made visible.
function SAPORMGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SAPORMGUI (see VARARGIN)

% Choose default command line output for SAPORMGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% This sets up the initial plot - only do when we are invisible
% so window can get raised using SAPORMGUI.
if strcmp(get(hObject,'Visible'),'off')
    % plot(rand(5));
end

% UIWAIT makes SAPORMGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SAPORMGUI_OutputFcn(hObject, eventdata, handles)
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
function RunSAP_Callback(hObject, eventdata, handles)


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
ORM = get(handles.ORM,'Value');
OffRM = get(handles.OffRM,'Value');
OnRM = get(handles.OnRM,'Value');
%---
MODS = [OnRM OffRM ORM];
%-------


%-------
doPlot = get(handles.doPlot,'Value');
PlotNum = str2double(get(handles.PlotNum,'String'));
doFluorPlot = get(handles.doFluorPlot,'Value');
FluorT = str2double(get(handles.FluorT,'String'));
%---
DOES = [doPlot PlotNum doFluorPlot FluorT];
%-------


%-------
doRev = get(handles.doRev,'Value');
Rev = str2double(get(handles.Rev,'String'));
%---
REVA = [doRev Rev];
%-------


%-------
doAMPAR = get(handles.doAMPAR,'Value');
G1RT = str2double(get(handles.G1RT,'String'));
G1ST = str2double(get(handles.G1ST,'String'));
AMPARrate = str2double(get(handles.AMPARrate,'String'));
AMPARnum = str2double(get(handles.AMPARnum,'String'));
%---
GLU = [doAMPAR G1RT AMPARrate AMPARnum G1ST];
%-------


%-------
GT1 = str2double(get(handles.GTon,'String'));
GT2 = str2double(get(handles.GToff,'String'));
GT3 = str2double(get(handles.LTPv,'String'));
gM = get(handles.gM,'Value');
%---
GT = [GT1 GT2 GT3 gM];
GTab = get(handles.GMASK,'Data');
%-------

%-------
doT1  = get(handles.do1h,'Value');
doT2  = get(handles.do4h,'Value');
doT3  = get(handles.do1d,'Value');
doT4  = get(handles.do1m,'Value');
doT5  = get(handles.do1y,'Value');
doT6  = str2double(get(handles.do1hT,'String'));
doT7  = str2double(get(handles.do4hT,'String'));
doT8  = str2double(get(handles.do1dT,'String'));
doT9  = str2double(get(handles.do1mT,'String'));
doT10 = str2double(get(handles.do1yT,'String'));
%---
doTs = [doT1 doT2 doT3 doT4 doT5 doT6 doT7 doT8 doT9 doT10];
%-------



% handles.output = SAPORM(LBR,TIME,SIZE,MODS,DOES,REVA,GLU,GT,GTab,doTs);
% handles.output = ORMATHFUN(LBR,TIME,SIZE,MODS,DOES,REVA,GLU,GT,GTab,doTs);
% handles.output = clustersim(LBR,TIME,SIZE,MODS,DOES,REVA,GLU,GT,GTab,doTs);
handles.output = ACTINFUN(LBR,TIME,SIZE,MODS,DOES,REVA,GLU,GT,GTab,doTs);

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








