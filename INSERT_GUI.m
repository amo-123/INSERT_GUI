function varargout = INSERT_GUI(varargin)
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @INSERT_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @INSERT_GUI_OutputFcn, ...
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


% --- Outputs from this function are returned to the command line.
function varargout = INSERT_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function INSERT_GUI_OpeningFcn(hObject, eventdata, handles, varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User actions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('_Convert');
% if ~isdir('_GenerateCorrections')
%     mkdir('_GenerateCorrections');
% end
addpath('_GenerateCorrections');


% Creating "Open" bushbuttons.
[Icon ColorMap] = imread('open.gif');
OpenIcon = ind2rgb(Icon,ColorMap);
set(handles.OpenWorkDir,'CData', OpenIcon);
set(handles.OpenXfile,'CData', OpenIcon);
set(handles.OpenYfile,'CData', OpenIcon);
set(handles.OpenFloodFile,'CData', OpenIcon);
set(handles.OpenPlayFile,'CData', OpenIcon);
set(handles.OpenTargetDir,'CData', OpenIcon);

% Initializing default values for main paremeters shown on GUI.
handles.NumHeads=10; % The same value to be set in Edit_NumHeads_CreateFcn too!!!!
handles.NumPM = 72;
handles.Head = 1;
handles.NumXlines = 19;
handles.NumYlines = 41;
handles.dx = 0.2;
handles.Dx = 2.222;

% Creating main path variables
% WorkDir is set to current folder, the rest is not set to any inital
% value.
handles.WorkDir = pwd;
set(handles.Edit_WorkDir,'String',handles.WorkDir,'TooltipString',handles.WorkDir);
handles.Xdir = '';
handles.Xfile = '';
handles.Ydir = '';
handles.Yfile = '';
handles.floodDir = '';
handles.floodFile = '';
handles.playDir = '';
handles.playFile = '';
handles.targetDir = '';

% By default, "Unif only" checkbox is not checked.
handles.unifOnly = 0;

% End of user actions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Choose default command line output for INSERT_GUI
handles.output = hObject;

guidata(hObject, handles);
OpenWorkDir_Callback(hObject, eventdata, handles)


% Do not delete this empty function!
function Figure_Main_CreateFcn(hObject, eventdata, handles)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Callback functions for GUI objects %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main parameters panel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Number of modules: NumHeads
function Edit_NumHeads_Callback(hObject, eventdata, handles)

num=str2double(get(hObject,'String'));

if (~isnan(num)) && (num>=1)
    handles.NumHeads = floor(num);

% If NumHeads has been changed, Number of Modules popup under calibration
% must be changed too.
strCell={};
for i=1:handles.NumHeads
    strCell{i} = num2str(i);
end
if get(handles.Popup_Head,'Value')>handles.NumHeads;
    set(handles.Popup_Head,'Value',1);
end
    set(handles.Popup_Head,'String',strCell);
else
    set(hObject,'String',handles.NumHeads);
end

guidata(hObject,handles);

function Edit_NumHeads_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

handles.NumHeads = 10; % Default vaue of NumHeads to be set here redundantly!!!!
guidata(hObject,handles);


% Number of PMs: NumPM
function Edit_NumPM_Callback(hObject, eventdata, handles)

num=str2double(get(hObject,'String'));

if (~isnan(num)) && (num>=1)
    handles.NumPM = floor(num);
else
    set(hObject,'String',handles.NumPM);
end

guidata(hObject,handles);


function Edit_NumPM_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% Pixel size: dx
function Edit_dx_Callback(hObject, eventdata, handles)

num=str2double(get(hObject,'String'));

if (~isnan(num)) && (num>=0)
    handles.dx = num;
else
    set(hObject,'String',handles.dx);
end

guidata(hObject,handles);

function Edit_dx_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% Setting main working directory
function OpenWorkDir_Callback(hObject, eventdata, handles)

dir=uigetdir(handles.WorkDir,'Select Root Working Directory...');
if (dir~=0)
    handles.WorkDir=dir;
end
set(handles.Edit_WorkDir,'String',handles.WorkDir,'TooltipString',handles.WorkDir);
guidata(hObject,handles);
ResetCalibDirs(hObject,handles);
handles=guidata(hObject);
ResetPlayDirs(hObject,handles);


function OpenWorkDir_CreateFcn(hObject, eventdata, handles)

function Edit_WorkDir_Callback(hObject, eventdata, handles)

function Edit_WorkDir_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calibration panel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Popup: selecting head
function Popup_Head_Callback(hObject, eventdata, handles)
str = get(hObject,'String');
val = get(hObject,'Value');
handles.Head = str2double(str(val));
guidata(hObject,handles);
ResetCalibDirs(hObject,handles);


function Popup_Head_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

strCell={};
for i=1:handles.NumHeads
    strCell{i} = num2str(i);
end
set(hObject,'String',strCell);


% NumXlines
function Edit_NumXlines_Callback(hObject, eventdata, handles)

num=str2double(get(hObject,'String'));

if (~isnan(num)) && (num>=1)
    handles.NumXlines = floor(num);
else
    set(hObject,'String',handles.NumXlines);
end

guidata(hObject,handles);

function Edit_NumXlines_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% NumYlines
function Edit_NumYlines_Callback(hObject, eventdata, handles)

num=str2double(get(hObject,'String'));

if (~isnan(num)) && (num>=1)
    handles.NumYlines = floor(num);
else
    set(hObject,'String',handles.NumYlines);
end

guidata(hObject,handles);

function Edit_NumYlines_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%Pitch: Dx
function Edit_Dx_Callback(hObject, eventdata, handles)

num=str2double(get(hObject,'String'));

if (~isnan(num)) && (num>=0)
    handles.Dx = num;
else
    set(hObject,'String',handles.Dx);
end

guidata(hObject,handles);

function Edit_Dx_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%Setting X phantom acquisition
function OpenXfile_Callback(hObject, eventdata, handles)
[file,dir]=uigetfile([handles.Xdir,'*.dat; *.data'],'Select X phantom acquisition...');
if (file~=0)
    handles.Xdir=dir;
    handles.Xfile = [dir,file];
end
set(handles.Edit_Xfile,'String',handles.Xfile,'TooltipString',handles.Xfile);
guidata(hObject,handles);

function Edit_Xfile_Callback(hObject, eventdata, handles)

function Edit_Xfile_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%Setting Y phantom acquisition
function OpenYfile_Callback(hObject, eventdata, handles)
[file,dir]=uigetfile([handles.Ydir,'*.dat; *.data'],'Select Y phantom acquisition...');
if (file~=0)
    handles.Ydir=dir;
    handles.Yfile = [dir,file];
end
set(handles.Edit_Yfile,'String',handles.Yfile,'TooltipString',handles.Yfile);
guidata(hObject,handles);

function Edit_Yfile_Callback(hObject, eventdata, handles)

function Edit_Yfile_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%Setting flood acquisition
function OpenFloodFile_Callback(hObject, eventdata, handles)
[file,dir]=uigetfile([handles.floodDir,'*.dat; *.data'],'Select flood acquisition...');
if (file~=0)
    handles.floodDir=dir;
    handles.floodFile = [dir,file];
end
set(handles.Edit_floodFile,'String',handles.floodFile,'TooltipString',handles.floodFile);
guidata(hObject,handles);


function Edit_floodFile_Callback(hObject, eventdata, handles)

function Edit_floodFile_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% Convert data pushbutton
function RunCalibConvert_Callback(hObject, eventdata, handles)
AbleButtons(handles,'off');


if handles.unifOnly == 0 % If "Unif only" is NOT checked
    
    try
        disp(' '); disp(' '); disp('Fnc_DataConverter >>>');
        [handles.Xdir,handles.Xfile] = ...
            Fnc_DataConverter(handles.Xfile,handles.WorkDir,handles.Head,'X');
        set(handles.Edit_Xfile,'String',handles.Xfile,'TooltipString',handles.Xfile);
        disp('X file converted.');
    catch
        disp('Could not convert X file.');
    end
    
    try
        disp(' '); disp(' '); disp('Fnc_DataConverter >>>');
        [handles.Ydir,handles.Yfile] = ...
            Fnc_DataConverter(handles.Yfile,handles.WorkDir,handles.Head,'Y');
        set(handles.Edit_Yfile,'String',handles.Yfile,'TooltipString',handles.Yfile);
        disp('Y file converted.');
    catch
        disp('Could not convert Y file.');
    end
    
end




try
    
    if strcmp(handles.floodFile((end-18):end), 'bulmaraw_ALL_F.data')
        
        disp(' '); disp(' '); disp('MexPlayDataConverter >>>');
        tic;
        MexPlayDataConverter(handles.floodFile,[handles.WorkDir,'\Corr\digitalraw'],handles.NumHeads);
        toc;
        ResetCalibDirs(hObject,handles);
        % Fnc_PlayDataConverter_04(handles.playFile,handles.targetDir,handles.NumHeads);
        disp('Flood file for all nodes is converted.');
    else
        disp(' '); disp(' '); disp('Fnc_DataConverter >>>');
        [handles.floodDir,handles.floodFile] = ...
            Fnc_DataConverter(handles.floodFile,handles.WorkDir,handles.Head,'F');
        set(handles.Edit_floodFile,'String',handles.floodFile,'TooltipString',handles.floodFile);
        disp('Flood file converted.');
    end
catch
    disp('Could not convert flood file.');
end

AbleButtons(handles,'on');
ResetCalibDirs(hObject,handles); % guidata(hObject,handles);


% Generate correction data pushbutton
function RunCalibration_Callback(hObject, eventdata, handles)
AbleButtons(handles,'off');

try
    % Kill the active figures
    r = groot;
    fig = r.Children;
    for i=1:length(fig)
        if strcmp(fig(i).Name,'INSERT_GUI') == 0
            delete(fig(i))
        end
    end
    
    if handles.unifOnly == 1  % If "Unif only" IS checked
        disp(' '); disp(' '); disp('Fnc_RefreshUnifCorr >>>');
        Fnc_RefreshUnifCorr(handles.Head, handles.NumXlines, handles.NumYlines, ...
                            handles.Dx, handles.dx, handles.NumPM, ...
                            handles.Xfile, handles.Yfile, handles.floodFile, ...
                            handles.WorkDir);
    else
        disp(' '); disp(' '); disp('Fnc_GenerateLRFandCorrections >>>');
        Fnc_GenerateLRFandCorrections(handles.Head, handles.NumXlines, handles.NumYlines, ...
                                      handles.Dx, handles.dx, handles.NumPM, ...
                                      handles.Xfile, handles.Yfile, handles.floodFile, ...
                                      handles.WorkDir);
    end
    
catch
    disp('Calibration failed.');
end

AbleButtons(handles,'on');


% Checkbox: Refresh Unif kalibration only
function Check_UnifOnly_Callback(hObject, eventdata, handles)

handles.unifOnly = get(hObject, 'value');

guidata(hObject,handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Play acquisition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Acquistion fie
function Edit_playFile_Callback(hObject, eventdata, handles)

function Edit_playFile_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function OpenPlayFile_Callback(hObject, eventdata, handles)
[file,dir]=uigetfile([handles.playFile,';.data'],'Select Acquisition...');
if (file~=0)
    handles.playDir=dir;
    handles.playFile = [dir,file];
end
set(handles.Edit_playFile,'String',handles.playFile,'TooltipString',handles.playFile);
guidata(hObject,handles);


% Target directory
function Edit_targetDir_Callback(hObject, eventdata, handles)

function Edit_targetDir_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function OpenTargetDir_Callback(hObject, eventdata, handles)
dir=uigetdir(handles.targetDir,'Select Target Directory...');
if (dir~=0)
    handles.targetDir=dir;
end
set(handles.Edit_targetDir,'String',handles.targetDir,'TooltipString',handles.targetDir);

guidata(hObject,handles);


% Run PlayConvert pushbutton
function RunPlayConvert_Callback(hObject, eventdata, handles)
AbleButtons(handles,'off');

try
    disp(' '); disp(' '); disp('MexPlayDataConverter >>>');
    tic;
    MexPlayDataConverter(handles.playFile,handles.targetDir,handles.NumHeads);
    % Fnc_PlayDataConverter_04(handles.playFile,handles.targetDir,handles.NumHeads);
    toc;
    disp('Play data converted.');
catch ME
    disp('Could not convert play data.');
    disp(ME.message);
end
AbleButtons(handles,'on');

% Run Play
function RunPlay_Callback(hObject, eventdata, handles)
AbleButtons(handles,'off');

try
    % Kill open figures
    r = groot;
    fig = r.Children;
    for i=1:length(fig)
        if strcmp(fig(i).Name,'INSERT_GUI') == 0
            delete(fig(i))
        end
    end

    disp(' '); disp(' '); disp('Fnc_PlayAcquisition >>>');
    Fnc_PlayAcquisition(handles);

catch ME
    disp('Could not play data.');
    disp(ME.message);
end


AbleButtons(handles,'on');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% General functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function ResetPlayDirs(hObject,handles)

handles.playDir = handles.WorkDir;
handles.playFile = [handles.WorkDir,'\*.*'];
set(handles.Edit_playFile,'String','[WorkDir]',...
   'TooltipString','Select acquisition!');


if ~exist([handles.WorkDir,'\Play'],'dir')
    
    AbleButtons(handles,'off')
    button = questdlg('Creating directory [WorkDir]\Play...','Confirmation','Yes','No','Yes');
    
    if strcmp(button,'Yes')
        mkdir([handles.WorkDir,'\Play']);
        handles.targetDir = [handles.WorkDir,'\Play'];
        set(handles.Edit_targetDir,'String','[WorkDir]\Play',...
           'TooltipString',handles.targetDir);
    else
        handles.targetDir = handles.WorkDir;
        set(handles.Edit_targetDir,'String','[WorkDir]',...
           'TooltipString',handles.targetDir);
    end
    AbleButtons(handles,'on')
else
      handles.targetDir = [handles.WorkDir,'\Play'];
      set(handles.Edit_targetDir,'String','[WorkDir]\Play',...
           'TooltipString',handles.targetDir);

end

guidata(hObject,handles);



function ResetCalibDirs(hObject,handles)


digitalrawDir = [handles.WorkDir,'\Corr\digitalraw\'];

digitalraw_X = ...
    [digitalrawDir,'digitalraw_H',num2str(handles.Head,'%02i'),'_X.dat'];
digitalraw_Y = ...
    [digitalrawDir,'digitalraw_H',num2str(handles.Head,'%02i'),'_Y.dat'];
digitalraw_F = ...
    [digitalrawDir,'digitalraw_H',num2str(handles.Head,'%02i'),'_F.dat'];
digitalraw = ...
    [digitalrawDir,'digitalraw_H',num2str(handles.Head,'%02i'),'.dat'];

bulmarawDir = [handles.WorkDir,'\Corr\bulmaraw\'];

bulmaraw_X = ...
    [bulmarawDir,'bulmaraw_H',num2str(handles.Head,'%02i'),'_X.data'];
bulmaraw_Y = ...
    [bulmarawDir,'bulmaraw_H',num2str(handles.Head,'%02i'),'_Y.data'];
bulmaraw_F = ...
    [bulmarawDir,'bulmaraw_H',num2str(handles.Head,'%02i'),'_F.data'];
bulmaraw_ALL_F = ...
    [bulmarawDir,'bulmaraw_ALL_F.data'];

if exist(digitalraw_X,'file')
    handles.Xdir = digitalrawDir;
    handles.Xfile = digitalraw_X;
    set(handles.Edit_Xfile,'String',['[WorkDir]',handles.Xfile((end-36):end)],...
        'TooltipString',handles.Xfile);
elseif exist(bulmaraw_X,'file')
    handles.Xdir = bulmarawDir;
    handles.Xfile = bulmaraw_X;
    set(handles.Edit_Xfile,'String',['[WorkDir]',handles.Xfile((end-33):end)],...
        'TooltipString',handles.Xfile);
else
    handles.Xdir = [handles.WorkDir,'\'];
    handles.Xfile = [handles.Xdir,'*.*'];
    set(handles.Edit_Xfile,'String','[WorkDir]',...
        'TooltipString',handles.Xfile);
end

if exist(digitalraw_Y,'file')
    handles.Ydir = digitalrawDir;
    handles.Yfile = digitalraw_Y;
    set(handles.Edit_Yfile,'String',['[WorkDir]',handles.Yfile((end-36):end)],...
        'TooltipString',handles.Yfile);
elseif exist(bulmaraw_Y,'file')
    handles.Ydir = bulmarawDir;
    handles.Yfile = bulmaraw_Y;
    set(handles.Edit_Yfile,'String',['[WorkDir]',handles.Yfile((end-33):end)],...
        'TooltipString',handles.Yfile);
else
    handles.Ydir = [handles.WorkDir,'\'];
    handles.Yfile = [handles.Ydir,'*.*'];
    set(handles.Edit_Yfile,'String','[WorkDir]',...
        'TooltipString',handles.Yfile);
end

if exist(digitalraw_F,'file')
    handles.floodDir = digitalrawDir;
    handles.floodFile = digitalraw_F;
    set(handles.Edit_floodFile,'String',['[WorkDir]',handles.floodFile((end-36):end)],...
        'TooltipString',handles.floodFile);
elseif exist(digitalraw,'file')
    handles.floodDir = digitalrawDir;
    handles.floodFile = digitalraw;
    set(handles.Edit_floodFile,'String',['[WorkDir]',handles.floodFile((end-34):end)],...
        'TooltipString',handles.floodFile);
elseif exist(bulmaraw_ALL_F,'file')
    handles.floodDir = bulmarawDir;
    handles.floodFile = bulmaraw_ALL_F;
    set(handles.Edit_floodFile,'String',['[WorkDir]',handles.floodFile((end-33):end)],...
        'TooltipString',handles.floodFile);
elseif exist(bulmaraw_F,'file')
    handles.floodDir = bulmarawDir;
    handles.floodFile = bulmaraw_F;
    set(handles.Edit_floodFile,'String',['[WorkDir]',handles.floodFile((end-33):end)],...
        'TooltipString',handles.floodFile);
else
    handles.floodDir = [handles.WorkDir,'\'];
    handles.floodFile = [handles.floodDir,'*.*'];
    set(handles.Edit_floodFile,'String','[WorkDir]',...
        'TooltipString',handles.floodFile);
end


guidata(hObject,handles);



% Enable/disable buttons while running
function AbleButtons(handles,onoff)

set(handles.Edit_NumHeads,'Enable',onoff);
set(handles.Edit_NumPM,'Enable',onoff);
set(handles.Edit_dx,'Enable',onoff);
set(handles.OpenWorkDir,'Enable',onoff);
set(handles.Popup_Head,'Enable',onoff);
set(handles.Edit_NumXlines,'Enable',onoff);
set(handles.Edit_NumYlines,'Enable',onoff);
set(handles.Edit_Dx,'Enable',onoff);
set(handles.OpenXfile,'Enable',onoff);
set(handles.OpenYfile,'Enable',onoff);
set(handles.OpenFloodFile,'Enable',onoff);
set(handles.RunCalibConvert,'Enable',onoff);
set(handles.RunCalibration,'Enable',onoff);
set(handles.OpenPlayFile,'Enable',onoff);
set(handles.OpenTargetDir,'Enable',onoff);
set(handles.RunPlayConvert,'Enable',onoff);
set(handles.RunPlay,'Enable',onoff);

pause(0.1);



% By pushing on the main window area, all buttons can be made enabled.
function Figure_Main_ButtonDownFcn(hObject, eventdata, handles)

AbleButtons(handles,'on');


% --- Executes on button press in pushbutton12. Batch file processing 
function batchproc(hObject, eventdata, handles)
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
folder =uigetdir([handles.playFile,';.data'],'Select Acquisition...');

files = dir(fullfile(folder, '*.data'));
batch = handles;
for i = 1:length(files)
    batch.playFile = [folder,'\',files(i).name];
    
    RunPlayConvert_Callback(handles.RunPlayConvert,[],batch);
    
    RunPlay_Callback(handles.RunPlay, [], batch);
    
    fid = load([handles.targetDir,'\corrected_images.mat']);
    save([handles.WorkDir,'\',num2str(files(i).name(1:end-5))],'fid');
    
    count = load([handles.targetDir,'\CountsAndEWC.mat']);
    save([handles.WorkDir,'\',num2str(files(i).name(1:end-5)),'_CountsAndEWC'],'count');
end

disp('Fin');
