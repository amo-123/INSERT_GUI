function varargout = corr_snakes(varargin)
% CORR_SNAKES MATLAB code for corr_snakes.fig
%      CORR_SNAKES, by itself, creates a new CORR_SNAKES or raises the existing
%      singleton*.
%
%      H = CORR_SNAKES returns the handle to a new CORR_SNAKES or the handle to
%      the existing singleton*.
%
%      CORR_SNAKES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CORR_SNAKES.M with the given input arguments.
%
%      CORR_SNAKES('Property','Value',...) creates a new CORR_SNAKES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before corr_snakes_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to corr_snakes_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help corr_snakes

% Last Modified by GUIDE v2.5 13-Mar-2019 11:32:05

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @corr_snakes_OpeningFcn, ...
                   'gui_OutputFcn',  @corr_snakes_OutputFcn, ...
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


% --- Executes just before corr_snakes is made visible.
function corr_snakes_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to corr_snakes (see VARARGIN)

% Choose default command line output for corr_snakes
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
init_cors();

% UIWAIT makes corr_snakes wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = corr_snakes_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in radiobutton1.                            [X]
function radiobutton1_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton1

cors = comdat('get','cors');
if get(hObject,'Value')
    cors.dim = 'X';
    set(handles.radiobutton2,'Value',0);
else
    cors.dim = 'Y';
    set(handles.radiobutton2,'Value',1);
end
comdat('set','cors',cors);

% --- Executes on button press in radiobutton2.                            [Y]
function radiobutton2_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton2

cors = comdat('get','cors');
if get(hObject,'Value')
    cors.dim = 'Y';
    set(handles.radiobutton1,'Value',0);
else
    cors.dim = 'X';
    set(handles.radiobutton1,'Value',1);
end
comdat('set','cors',cors);

% --- Executes on slider movement.                                         [X]
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

cors = comdat('get','cors');
cors.ix_snk = round( get(hObject,'Value') );
set(handles.edit1,'String',num2str(cors.ix_snk))
comdat('set','cors',cors);

% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.                                         [Y]
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

cors = comdat('get','cors');
cors.iy_snk = round( get(hObject,'Value') );
set(handles.edit2,'String',num2str(cors.iy_snk))
comdat('set','cors',cors);

% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit1_Callback(hObject, eventdata, handles)                      %[X]
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


function edit2_Callback(hObject, eventdata, handles)                      %[Y]
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


% --- Executes on button press in pushbutton1.                             [Finish]
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cors = comdat('get','cors');

delete(cors.n_fig);
%close(corr_snakes);



% --- Executes on button press in pushbutton2.                             [Go]
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

cors = comdat('get','cors');
LineData = comdat('get','LineData');
ImageData = comdat('get','ImageData');
if cors.dim == 'X'
    i_snk = cors.ix_snk;
else 
    i_snk = cors.iy_snk;
end
LineMin = find(LineData(i_snk,:,2)>0,1,'first');
LineMax = find(LineData(i_snk,:,2)>0,1,'last');
Axis= LineMin:LineMax;
figure(cors.n_fig);
imagesc(-ImageData), hold on
switch cors.dim
    case 'X'
        plot(LineData(i_snk,Axis,2),Axis,'w')
    case 'Y'
        plot(Axis,LineData(i_snk,Axis,2),'w')
end
hold off;
zoom(1.8);

[x, y] = ginput;

NewLine = LineData(i_snk,:,2);

try
    switch cors.dim
        case 'X'
            jY = round(min(y)):round(max(y));
            NewLine(jY) = interp1(y,x,jY);
        case 'Y'
            jX = round(min(x)):round(max(x));
            NewLine(jX) = interp1(x,y,jX);
    end
    
    
    figure(cors.n_fig);
    imagesc(-ImageData), hold on
    switch cors.dim
        case 'X'
            plot(NewLine(Axis),Axis,'w')
        case 'Y'
            plot(Axis,NewLine(Axis),'w')
    end
    hold off;
    zoom(1.8);
    
    cors.NewSnake = NewLine;
    cors.i_new = i_snk;
    
    comdat('set','cors',cors);
catch
    errordlg('Please pick more data ponits')
end





%===
function init_cors()

fig=figure; 

cors.n_fig = fig.Number;
cors.dim = 'X';
cors.ix_snk = 1;
cors.iy_snk = 1;
cors.NewSnake = [];

comdat('set','cors',cors)


% --- Executes on button press in pushbutton3.                             Confirm
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

cors = comdat('get','cors');
LineData = comdat('get','LineData');

LineData(cors.i_new,:,2) = cors.NewSnake;

comdat('set', 'LineData', LineData);


% --- Executes on button press in pushbutton4.                             View
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

cors = comdat('get','cors');
LineData = comdat('get','LineData');
ImageData = comdat('get','ImageData');
if cors.dim == 'X'
    n_snk = 19;
else 
    n_snk = 41;
end

figure(cors.n_fig);
imagesc(-ImageData), hold on
for i = 1:n_snk
    LineMin = find(LineData(i,:,2)>0,1,'first');
    LineMax= find(LineData(i,:,2)>0,1,'last');
    Axis= LineMin:LineMax;
    switch cors.dim
        case 'X'
            plot(LineData(i,Axis,2),Axis,'w')
        case 'Y'
            plot(Axis,LineData(i,Axis,2),'w')
    end
end
hold off;
zoom(1.8);


% --- Executes on button press in pushbutton5.                             Help
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
msgbox({'Select the axis and snake number to edit';...
    'Press Go and select points to draw a new path for the chosen snake';...
    'When finished selecting points, press enter key';...
    'New snake will be shown, confirm to save new snake';...
    'View all current snakes before selecting finished'});


% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in pushbutton6.                             Cut Off 
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cors = comdat('get','cors');
LineData= comdat('get','LineData');
ImageData = comdat('get','ImageData');

SolidSnk = LineData(:,:,2);

if cors.dim == 'X'
    n_snk = 19;
else 
    n_snk = 41;
end
mid = zeros(n_snk,1);
for i = 1:n_snk
    LineMin = find(SolidSnk(i,:)>0,1,'first');
    LineMax= find(SolidSnk(i,:)>0,1,'last');
    Axis= LineMin:LineMax;
    mid(i) = median(Axis);
end

TrueMid = ceil(mean(mid));
figure(cors.n_fig);
imagesc(-ImageData), hold on
for i = 1:n_snk  
    LineMin = find(SolidSnk(i,:)>0,1,'first');
    LineMax= find(SolidSnk(i,:)>0,1,'last');
    switch cors.dim
        case 'X'
            NewMin = TrueMid - 250;
            NewMax = TrueMid + 250;
            NewAxis = NewMin:NewMax;
            MGS = 1:length(SolidSnk(i,:));
            ind = MGS < NewMin;
            SolidSnk(i,ind) = 0;
            ind = MGS > NewMax;
            SolidSnk(i,ind) = 0;
            try
                ind = SolidSnk(i,:)> NewMin && SolidSnk(i,:) < LineMin;
                SolidSnk(i,ind) = SolidSnk(i,LineMin);
                ind = SolidSnk(i,:)> LineMax && SolidSnk(i,:) < NewMax;
                SolidSnk(i,ind) = SolidSnk(i,LineMin);
            catch
                
            end
            plot(SolidSnk(i,NewAxis),NewAxis,'.w')
        case 'Y'
            NewMin = TrueMid - 120;
            NewMax = TrueMid + 120;
            NewAxis = NewMin:NewMax;
            MGS = 1:length(SolidSnk(i,:));
            ind = MGS < NewMin;
            SolidSnk(i,ind) = 0;
            ind = MGS > NewMax;
            SolidSnk(i,ind) = 0;
            try
                ind = MGS < LineMin & MGS > NewMin;
                SolidSnk(i,ind) = SolidSnk(i,LineMin);
                ind = MGS > LineMax & MGS < NewMax;
                SolidSnk(i,ind) = SolidSnk(i,LineMax);
            catch
                
            end
            plot(NewAxis, SolidSnk(i,NewAxis),'.w')            
    end
end
hold off;
zoom(1.8);


LineData(:,:,2) = SolidSnk;

comdat('set', 'LineData', LineData);


