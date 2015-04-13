function varargout = ajusta_impedancias(varargin)
% AJUSTA_IMPEDANCIAS MATLAB code for ajusta_impedancias.fig
%      AJUSTA_IMPEDANCIAS, by itself, creates a new AJUSTA_IMPEDANCIAS or raises the existing
%      singleton*.
%
%      H = AJUSTA_IMPEDANCIAS returns the handle to a new AJUSTA_IMPEDANCIAS or the handle to
%      the existing singleton*.
%
%      AJUSTA_IMPEDANCIAS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in AJUSTA_IMPEDANCIAS.M with the given input arguments.
%
%      AJUSTA_IMPEDANCIAS('Property','Value',...) creates a new AJUSTA_IMPEDANCIAS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ajusta_impedancias_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ajusta_impedancias_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ajusta_impedancias

% Last Modified by GUIDE v2.5 25-Apr-2013 10:51:43

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ajusta_impedancias_OpeningFcn, ...
                   'gui_OutputFcn',  @ajusta_impedancias_OutputFcn, ...
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


% --- Executes just before ajusta_impedancias is made visible.
function ajusta_impedancias_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ajusta_impedancias (see VARARGIN)
handles.var.plano = 1;
handles.var.codigo = 1;
handles.var.res_list = get(handles.popup_sel_res,'String');
handles.var.active_res = 1;
handles.var.exp_max_iter = 3;

cla(handles.ReZ);
box(handles.ReZ,'on');
hold(handles.ReZ,'all');
% Create xlabel
xlabel(handles.ReZ, 'f [GHz]');
% Create ylabel
ylabel(handles.ReZ, 'Re(Z) [SI]');

cla(handles.ImZ);
box(handles.ImZ,'on');
hold(handles.ImZ,'all');
% Create xlabel
xlabel(handles.ImZ,'f [GHz]');
% Create ylabel
ylabel(handles.ImZ,'Im(Z) [SI]');

num_res = 8;
handles.var.num_res = num_res;
handles.var.Rs = zeros(1,num_res);
handles.var.wr = 15*2*pi*1e9*ones(1,num_res);
handles.var.Q  = 100*ones(1,num_res);
handles.var.Rs_min = zeros(1,num_res);
handles.var.wr_min = 1*2*pi*1e9*ones(1,num_res);
handles.var.Q_min  = 1*ones(1,num_res);
handles.var.Rs_max = 100*ones(1,num_res);
handles.var.wr_max = 100*2*pi*1e9*ones(1,num_res);
handles.var.Q_max  = 1000*ones(1,num_res);
handles.var.Zfit = [];
handles.var.Z  = [];
handles.var.w  = [];
% Define a semente dos numeros aleatorios:
RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)));

handles.var.path = pwd;

% Choose default command line output for ajusta_impedancias
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ajusta_impedancias wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ajusta_impedancias_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%% CARREGA ARQUIVO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on selection change in popup_plano.
function popup_plano_Callback(hObject, eventdata, handles)
% hObject    handle to popup_plano (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes during object creation, after setting all properties.
function popup_plano_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_plano (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
string = {'Longitudinal', 'Vertical', 'Horizontal'};
set(hObject,'String',string);
set(hObject,'Value',1);



% --- Executes on selection change in popup_codigo.
function popup_codigo_Callback(hObject, eventdata, handles)
% hObject    handle to popup_codigo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes during object creation, after setting all properties.
function popup_codigo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_codigo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
string = {'GdfidL', 'CST', 'ACE3P', 'ECHO'};
set(hObject,'String',string);
set(hObject,'Value',1);


% --- Executes on button press in pushbutton_diretorio.
function pushbutton_diretorio_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_diretorio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% load impedance
path = uigetdir('.','Em qual Diretorio estao os arquivos?');
set(handles.text_diretorio,'String',path);
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function text_diretorio_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_diretorio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(hObject,'String','Aqui aparecera o caminho do diretorio escolhido');


% --- Executes on button press in pushbutton_carrega_imped.
function pushbutton_carrega_imped_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_carrega_imped (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

plano = get(handles.popup_plano,'String');
pl = get(handles.popup_plano,'Value');
plano = plano{pl};

codigo = get(handles.popup_codigo,'String');
cdg = get(handles.popup_codigo,'Value');
codigo = codigo{cdg};

path = get(handles.text_diretorio,'String');

[Z w] = lnls_load_impedance(path, plano, codigo);
handles.var.Z = Z;
handles.var.w = w;

handles.var.plano  = pl;
handles.var.codigo = cdg;
handles.var.path   = path;

% Muda de ressonador, dependendo do plano escolhido
if handles.var.plano ==1
    handles.var.function = @lnls_calc_impedance_longitudinal_resonator;
else
    handles.var.function = @lnls_calc_impedance_transverse_resonator;
end

% Uncomment the following line to preserve the X-limits of the axes
xlim(handles.ReZ,[0 w(end)]/(1e9*2*pi));
cla(handles.ReZ);
% Create plot
handles.plotReZ = plot(w(floor(end/2):end)/(1e9*2*pi),real(Z(floor(end/2):end)),'Parent',handles.ReZ,'Color','b');

% Uncomment the following line to preserve the X-limits of the axes
xlim(handles.ImZ,[0 w(end)]/(1e9*2*pi));
cla(handles.ImZ);
% Create plot
handles.plotImZ = plot(w(floor(end/2):end)/(1e9*2*pi),imag(Z(floor(end/2):end)),'Parent',handles.ImZ,'Color','b');

guidata(hObject, handles);


%% NUMERO DE RESSONADORES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function edit_num_res_Callback(hObject, eventdata, handles)
% hObject    handle to edit_num_res (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_num_res as text
%        str2double(get(hObject,'String')) returns contents of edit_num_res as a double
num_res = str2double(get(hObject,'String'));

if num_res < 1
    set(hObject,'String','Invalid');
    return;
end
handles.var.num_res = num_res;

Rs = handles.var.Rs;
wr = handles.var.wr;
Q  = handles.var.Q;
Rs_max = handles.var.Rs_max;
wr_max = handles.var.wr_max;
Q_max  = handles.var.Q_max;
Rs_min = handles.var.Rs_min;
wr_min = handles.var.wr_min;
Q_min  = handles.var.Q_min;
size_antigo = length(Rs);

if size_antigo > handles.var.num_res;
    handles.var.Rs = Rs(1:num_res);
    handles.var.wr = wr(1:num_res);
    handles.var.Q  = Q(1:num_res);
    handles.var.Rs_max = Rs_max(1:num_res);
    handles.var.wr_max = wr_max(1:num_res);
    handles.var.Q_max  = Q_max(1:num_res);
    handles.var.Rs_min = Rs_min(1:num_res);
    handles.var.wr_min = wr_min(1:num_res);
    handles.var.Q_min  = Q_min(1:num_res);
else
    handles.var.Rs = [Rs zeros(1,num_res-size_antigo)];
    handles.var.wr = [wr 50e9*ones(1,num_res-size_antigo)];
    handles.var.Q  = [Q  100*ones(1,num_res-size_antigo)];
    handles.var.Rs_max = [Rs_max 100*ones(1,num_res-size_antigo)];
    handles.var.wr_max = [wr_max 2*pi*100e9*ones(1,num_res-size_antigo)];
    handles.var.Q_max  = [Q_max  1000*ones(1,num_res-size_antigo)];
    handles.var.Rs_min = [Rs_min zeros(1,num_res-size_antigo)];
    handles.var.wr_min = [wr_min 2*pi*1e9*ones(1,num_res-size_antigo)];
    handles.var.Q_min  = [Q_min  1*ones(1,num_res-size_antigo)];
end

for m=1:handles.var.num_res
    if m==1
        handles.var.res_list = {'res 1'};
    else
        handles.var.res_list = [handles.var.res_list {sprintf('res %d',m)}];
    end
end

set(handles.popup_sel_res,'String', handles.var.res_list);
res_act = 1;
handles.var.active_res = res_act;
set(handles.popup_sel_res,'Value', res_act);
set(handles.edit_fr_min,'String',num2str(handles.var.wr_min(res_act)/(1e9*2*pi)));
set(handles.edit_Rs_min,'String',num2str(handles.var.Rs_min(res_act)));
set(handles.edit_Q_min,'String',num2str(handles.var.Q_min(res_act)));
set(handles.edit_fr_max,'String',num2str(handles.var.wr_max(res_act)/(1e9*2*pi)));
set(handles.edit_Rs_max,'String',num2str(handles.var.Rs_max(res_act)));
set(handles.edit_Q_max,'String',num2str(handles.var.Q_max(res_act)));
set(handles.edit_fr,'String',num2str(handles.var.wr(res_act)/(1e9*2*pi)));
set(handles.edit_Rs,'String',num2str(handles.var.Rs(res_act)));
set(handles.edit_Q,'String',num2str(handles.var.Q(res_act)));
set(handles.slider_fr,'Max',handles.var.wr_max(res_act)/(1e9*2*pi));
set(handles.slider_Rs,'Max',handles.var.Rs_max(res_act));
set(handles.slider_Q,'Max',handles.var.Q_max(res_act));
set(handles.slider_fr,'Value',handles.var.wr(res_act)/(1e9*2*pi));
set(handles.slider_Rs,'Value',handles.var.Rs(res_act));
set(handles.slider_Q,'Value',handles.var.Q(res_act));
set(handles.slider_fr,'Min',handles.var.wr_min(res_act)/(1e9*2*pi));
set(handles.slider_Rs,'Min',handles.var.Rs_min(res_act));
set(handles.slider_Q,'Min',handles.var.Q_min(res_act));

handles = calcula_plota_impedancia(handles);

guidata(hObject, handles);

function handles = calcula_plota_impedancia(handles)

handles.var.Zfit = handles.var.function(handles.var.Rs, handles.var.Q, handles.var.wr, handles.var.w);
vec_residue = [real(handles.var.Z - handles.var.Zfit) imag(handles.var.Z - handles.var.Zfit)]';
handles.var.residue = vec_residue'*vec_residue;
set(handles.text_residuo_atual,'String',sprintf('%6.4g',handles.var.residue));

get(handles.ReZ,'Children');
try
    delete(handles.plotReZfit);
    delete(handles.plotImZfit);
catch
end
handles.plotReZfit = plot(handles.var.w(floor(end/2):end)/(1e9*2*pi),real(handles.var.Zfit(floor(end/2):end)),'Parent',handles.ReZ,'Color','r');
handles.plotImZfit = plot(handles.var.w(floor(end/2):end)/(1e9*2*pi),imag(handles.var.Zfit(floor(end/2):end)),'Parent',handles.ImZ,'Color','r');
drawnow;

% --- Executes during object creation, after setting all properties.
function edit_num_res_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_num_res (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

set(hObject,'String',num2str(8));
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



%% AJUSTE INDIVIDUAL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on selection change in popup_sel_res.
function popup_sel_res_Callback(hObject, eventdata, handles)
% hObject    handle to popup_sel_res (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_sel_res contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_sel_res
handles.var.active_res = get(hObject,'Value');
res_act = handles.var.active_res;
set(handles.edit_fr_min,'String',num2str(handles.var.wr_min(res_act)/(1e9*2*pi)));
set(handles.edit_Rs_min,'String',num2str(handles.var.Rs_min(res_act)));
set(handles.edit_Q_min,'String',num2str(handles.var.Q_min(res_act)));
set(handles.edit_fr_max,'String',num2str(handles.var.wr_max(res_act)/(1e9*2*pi)));
set(handles.edit_Rs_max,'String',num2str(handles.var.Rs_max(res_act)));
set(handles.edit_Q_max,'String',num2str(handles.var.Q_max(res_act)));
set(handles.edit_fr,'String',num2str(handles.var.wr(res_act)/(1e9*2*pi)));
set(handles.edit_Rs,'String',num2str(handles.var.Rs(res_act)));
set(handles.edit_Q,'String',num2str(handles.var.Q(res_act)));
set(handles.slider_fr,'Max',handles.var.wr_max(res_act)/(1e9*2*pi));
set(handles.slider_Rs,'Max',handles.var.Rs_max(res_act));
set(handles.slider_Q,'Max',handles.var.Q_max(res_act));
set(handles.slider_fr,'Value',handles.var.wr(res_act)/(1e9*2*pi));
set(handles.slider_Rs,'Value',handles.var.Rs(res_act));
set(handles.slider_Q,'Value',handles.var.Q(res_act));
set(handles.slider_fr,'Min',handles.var.wr_min(res_act)/(1e9*2*pi));
set(handles.slider_Rs,'Min',handles.var.Rs_min(res_act));
set(handles.slider_Q,'Min',handles.var.Q_min(res_act));
guidata(hObject,handles);
% --- Executes during object creation, after setting all properties.
function popup_sel_res_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_sel_res (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
for m=1:8
    if m==1
        string = {'res 1'};
    else
        string = [string {sprintf('res %d',m)}];
    end
end
set(hObject,'String', string);



function edit_fr_max_Callback(hObject, eventdata, handles)
% hObject    handle to edit_fr_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.var.wr_max(handles.var.active_res) = (1e9*2*pi)*str2double(get(hObject,'String'));
set(handles.slider_fr,'Max', handles.var.wr_max(handles.var.active_res)/(1e9*2*pi));
guidata(hObject, handles);
% --- Executes during object creation, after setting all properties.
function edit_fr_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_fr_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(hObject,'String',num2str(100));
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_fr_min_Callback(hObject, eventdata, handles)
% hObject    handle to edit_fr_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.var.wr_min(handles.var.active_res) = (1e9*2*pi)*str2double(get(hObject,'String'));
set(handles.slider_fr,'Min', handles.var.wr_min(handles.var.active_res)/(1e9*2*pi));
guidata(hObject, handles);
% --- Executes during object creation, after setting all properties.
function edit_fr_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_fr_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(hObject,'String',num2str(1));
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_fr_Callback(hObject, eventdata, handles)
% hObject    handle to edit_fr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.var.wr(handles.var.active_res) = (1e9*2*pi)*str2double(get(hObject,'String'));
set(handles.slider_fr,'Value', handles.var.wr(handles.var.active_res)/(1e9*2*pi));
handles = calcula_plota_impedancia(handles);
guidata(hObject, handles);
% --- Executes during object creation, after setting all properties.
function edit_fr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_fr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(hObject,'String',num2str(15));
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Rs_max_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Rs_max (see GCBO)
% eventdata  reserved - to be defined in a future veion of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.var.Rs_max(handles.var.active_res) = str2double(get(hObject,'String'));
set(handles.slider_Rs,'Max', handles.var.Rs_max(handles.var.active_res));
guidata(hObject, handles);
% --- Executes during object creation, after setting all properties.
function edit_Rs_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Rs_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(hObject,'String',num2str(100));
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_Rs_min_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Rs_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.var.Rs_min(handles.var.active_res) = str2double(get(hObject,'String'));
set(handles.slider_Rs,'Min', handles.var.Rs_min(handles.var.active_res));
guidata(hObject, handles);
% --- Executes during object creation, after setting all properties.
function edit_Rs_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Rs_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(hObject,'String',num2str(0));
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_Rs_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Rs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.var.Rs(handles.var.active_res) = str2double(get(hObject,'String'));
set(handles.slider_Rs,'Value', handles.var.Rs(handles.var.active_res));
handles = calcula_plota_impedancia(handles);
guidata(hObject, handles);
% --- Executes during object creation, after setting all properties.
function edit_Rs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Rs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(hObject,'String',num2str(0));
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Q_max_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Q_max (see GCBO)
% eventdata  reserved - to be defined in a future veion of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.var.Q_max(handles.var.active_res) = str2double(get(hObject,'String'));
set(handles.slider_Q,'Max', handles.var.Q_max(handles.var.active_res));
guidata(hObject, handles);
% --- Executes during object creation, after setting all properties.
function edit_Q_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Q_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(hObject,'String',num2str(1000));
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_Q_min_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Q_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.var.Q_min(handles.var.active_res) = str2double(get(hObject,'String'));
set(handles.slider_Q,'Min', handles.var.Q_min(handles.var.active_res));
guidata(hObject, handles);
% --- Executes during object creation, after setting all properties.
function edit_Q_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Q_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(hObject,'String',num2str(1));
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_Q_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Q (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.var.Q(handles.var.active_res) = str2double(get(hObject,'String'));
set(handles.slider_Q,'Value', handles.var.Q(handles.var.active_res));
handles = calcula_plota_impedancia(handles);
guidata(hObject, handles);
% --- Executes during object creation, after setting all properties.
function edit_Q_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Q (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(hObject,'String',num2str(100));
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on slider movement.
function slider_fr_Callback(hObject, eventdata, handles)
% hObject    handle to slider_fr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
val = handles.var.active_res;
handles.var.wr(val) = (1e9*2*pi)*get(hObject,'Value');
set(handles.edit_fr,'String',num2str(handles.var.wr(val)/(1e9*2*pi)));
handles = calcula_plota_impedancia(handles);
guidata(hObject, handles);
% --- Executes during object creation, after setting all properties.
function slider_fr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_fr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(hObject,'Max',100);
set(hObject,'Value',15);
set(hObject,'Min',1);
% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



% --- Executes on slider movement.
function slider_Rs_Callback(hObject, eventdata, handles)
% hObject    handle to slider_Rs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
val = handles.var.active_res;
handles.var.Rs(val) = get(hObject,'Value');
set(handles.edit_Rs,'String',num2str(handles.var.Rs(val)));
handles = calcula_plota_impedancia(handles);
guidata(hObject, handles);
% --- Executes during object creation, after setting all properties.
function slider_Rs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_Rs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(hObject,'Max',100);
set(hObject,'Value',0);
set(hObject,'Min',0);
% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



% --- Executes on slider movement.
function slider_Q_Callback(hObject, eventdata, handles)
% hObject    handle to slider_Q (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
val = handles.var.active_res;
handles.var.Q(val) = get(hObject,'Value');
set(handles.edit_Q,'String',num2str(handles.var.Q(val)));
handles = calcula_plota_impedancia(handles);
guidata(hObject, handles);
% --- Executes during object creation, after setting all properties.
function slider_Q_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_Q (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(hObject,'Max',1000);
set(hObject,'Value',100);
set(hObject,'Min',1);
% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



% --- Executes on button press in pushbutton_Otimiza_Indiv.
function pushbutton_Otimiza_Indiv_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Otimiza_Indiv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(hObject,'Enable','off');
string = get(hObject,'String');
set(hObject,'String','Otimizando...');
drawnow;

Rs = handles.var.Rs;
wr = handles.var.wr;
Q  = handles.var.Q;
w  = handles.var.w;
num_max_iter = 10^handles.var.exp_max_iter;
residue = handles.var.residue;
res_act = handles.var.active_res;

par = [Rs(res_act) Q(res_act) wr(res_act)];

assignin('base','estado',[Rs;wr;Q]);

% maximo erro relativo permitido:
err = 0.01;
for i=1:num_max_iter
    error = (err*2*(rand(size(par))-0.5));
    new_par = par.*(1+error);
    Rs(res_act) = new_par(1);
    Q(res_act)  = new_par(2);
    wr(res_act) = new_par(3);
    new_Z  = handles.var.function(Rs, Q, wr, w);
    vec_new_residue = [real(handles.var.Z - new_Z) imag(handles.var.Z - new_Z)]';
    new_residue = vec_new_residue'*vec_new_residue;
    change = residue - new_residue;
    if change > 0
        residue = new_residue;
        par = new_par;
        handles.var.Rs = Rs;
        handles.var.wr = wr;
        handles.var.Q  = Q;
        handles = calcula_plota_impedancia(handles);
    end
end
set(handles.edit_fr,'String',num2str(handles.var.wr(res_act)/(1e9*2*pi)));
set(handles.edit_Rs,'String',num2str(handles.var.Rs(res_act)));
set(handles.edit_Q,'String',num2str(handles.var.Q(res_act)));
set(handles.slider_fr,'Value',handles.var.wr(res_act)/(1e9*2*pi));
set(handles.slider_Rs,'Value',handles.var.Rs(res_act));
set(handles.slider_Q,'Value',handles.var.Q(res_act));

set(hObject,'Enable','on');
set(hObject,'String',string);

guidata(hObject,handles);


% --- Executes on button press in push_recupera_estado.
function push_recupera_estado_Callback(hObject, eventdata, handles)
% hObject    handle to push_recupera_estado (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
estado = evalin('base','estado');
handles.var.Rs = estado(1,:);
handles.var.wr = estado(2,:);
handles.var.Q  = estado(3,:);

handles = calcula_plota_impedancia(handles);
res_act = handles.var.active_res;
set(handles.edit_fr,'String',num2str(handles.var.wr(res_act)/(1e9*2*pi)));
set(handles.edit_Rs,'String',num2str(handles.var.Rs(res_act)));
set(handles.edit_Q,'String',num2str(handles.var.Q(res_act)));
set(handles.slider_fr,'Value',handles.var.wr(res_act)/(1e9*2*pi));
set(handles.slider_Rs,'Value',handles.var.Rs(res_act));
set(handles.slider_Q,'Value',handles.var.Q(res_act));
guidata(hObject,handles);


function edit_max_itera_indiv_Callback(hObject, eventdata, handles)
% hObject    handle to edit_max_itera_indiv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.var.exp_max_iter = str2double(get(hObject,'String'));
guidata(hObject,handles);
% --- Executes during object creation, after setting all properties.
function edit_max_itera_indiv_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_max_itera_indiv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(hObject,'String', num2str(3));
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% OTIMIZADOR GERAL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_target_residue_Callback(hObject, eventdata, handles)
% hObject    handle to edit_target_residue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function edit_target_residue_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_target_residue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(hObject,'String','100');
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_max_iter_Callback(hObject, eventdata, handles)
% hObject    handle to edit_max_iter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit_max_iter as text
%        str2double(get(hObject,'String')) returns contents of edit_max_iter as a double

% --- Executes during object creation, after setting all properties.
function edit_max_iter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_max_iter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(hObject,'String','1');
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in push_comeca_svd.
function push_comeca_svd_Callback(hObject, eventdata, handles)
% hObject    handle to push_comeca_svd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% max_it = str2double(get(handles.edit_max_iter,'String'));
% Z  = handles.var.Z;
% w  = handles.var.w;
% Zfit = handles.var.Zfit;
% 
% target_residue = str2double(get(handles.edit_target_residue,'String'));
% vec_control = [ handles.var.Rs handles.var.Q handles.var.wr]';
% vec_residue = [real(Z-Zfit) imag(Z-Zfit)]';
% nr_sing_values = length(vec_control);
% 
% for i=1:max_it
%     M = calcula_matriz(vec_control,w,Z,handles);
%     [U,S,V] = svd(M,'econ');
%     iS = diag(1./diag(S));
%     diS = diag(iS);
%     diS(nr_sing_values+1:end) = 0;
%     iS = diag(diS);
%     CM = (V*iS*U');
%     delta_vec_control = CM * vec_residue;
%     vec_control = vec_control + delta_vec_control;
%     handles.var.Rs = vec_control(1:end/3)';
%     handles.var.Q  = vec_control(end/3+1:2*end/3)';
%     handles.var.wr = vec_control(2*end/3+1:end)';
%     handles = calcula_plota_impedancia(handles);
%     guidata(hObject, handles);
%     if handles.var.residue < target_residue
%         break;
%     end
% end
% function M = calcula_matriz(vec_control,w,Z,handles)
% M = zeros(length([real(Z) imag(Z)]), length(vec_control));
% new_vec_control = vec_control;
% for i=1:length(vec_control)
%     new_vec_control(i) = vec_control(i)*(1 - 0.02);
%     Rs = new_vec_control(1:end/3)';
%     Q  = new_vec_control(end/3+1:2*end/3)';
%     wr = new_vec_control(2*end/3+1:end)';
%     Zfit1= handles.var.function(Rs, Q, wr,w);
%     vec_target1 = [real(Zfit1) imag(Zfit1)]';
%     new_vec_control(i) = vec_control(i)*(1 + 0.04);
%     Rs = new_vec_control(1:end/3)';
%     Q  = new_vec_control(end/3+1:2*end/3)';
%     wr = new_vec_control(2*end/3+1:end)';
%     Zfit2= handles.var.function(Rs, Q, wr,w);
%     vec_target2 = [real(Zfit2) imag(Zfit2)]';
%     new_vec_control(i) = vec_control(i)*(1 - 0.02);
%     M(:,i) = (vec_target2 - vec_target1) / (vec_control(i)*0.04);
% end
set(hObject,'Enable','off');
string = get(hObject,'String');
set(hObject,'String','Otimizando...');
drawnow;

Rs = handles.var.Rs;
wr = handles.var.wr;
Q  = handles.var.Q;
w  = handles.var.w;
num_max_iter = 10^handles.var.exp_max_iter;
residue = handles.var.residue;
res_act = handles.var.active_res;



assignin('base','estado',[Rs;wr;Q]);
% maximo erro relativo permitido:
err = 0.01;
for jj=1:length(Rs)
    Rs = handles.var.Rs;
    wr = handles.var.wr;
    Q  = handles.var.Q;
    par = [Rs(jj) Q(jj) wr(jj)];
    for i=1:num_max_iter
        error = (err*2*(rand(size(par))-0.5));
        new_par = par.*(1+error);
        Rs(jj) = new_par(1);
        Q(jj)  = new_par(2);
        wr(jj) = new_par(3);
        new_Z  = handles.var.function(Rs, Q, wr, w);
        vec_new_residue = [real(handles.var.Z - new_Z) imag(handles.var.Z - new_Z)]';
        new_residue = vec_new_residue'*vec_new_residue;
        change = residue - new_residue;
        if change > 0
            residue = new_residue;
            par = new_par;
            handles.var.Rs = Rs;
            handles.var.wr = wr;
            handles.var.Q  = Q;
            handles = calcula_plota_impedancia(handles);
        end
    end
end
set(handles.edit_fr,'String',num2str(handles.var.wr(res_act)/(1e9*2*pi)));
set(handles.edit_Rs,'String',num2str(handles.var.Rs(res_act)));
set(handles.edit_Q,'String',num2str(handles.var.Q(res_act)));
set(handles.slider_fr,'Value',handles.var.wr(res_act)/(1e9*2*pi));
set(handles.slider_Rs,'Value',handles.var.Rs(res_act));
set(handles.slider_Q,'Value',handles.var.Q(res_act));

set(hObject,'Enable','on');
set(hObject,'String',string);

guidata(hObject,handles);


%% CARREGA AJUSTE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in push_carrega_arquivo.
function push_carrega_arquivo_Callback(hObject, eventdata, handles)
% hObject    handle to push_carrega_arquivo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
string = get(handles.popup_plano,'String');
string2 = get(handles.popup_codigo,'String');
fileName = ['Zfit_' string{handles.var.plano} '_' string2{handles.var.codigo} '.dat'];
fileName = fullfile(handles.var.path,fileName);
A = importdata(fileName)';
handles.var.Rs = A.data(:,1)';
handles.var.wr = A.data(:,2)'*1e9;
handles.var.Q  = A.data(:,3)';
handles.var.num_res = length(handles.var.Rs);

set(handles.edit_num_res,'String',num2str(handles.var.num_res));

for m=1:handles.var.num_res
    if m==1
        handles.var.res_list = {'res 1'};
    else
        handles.var.res_list = [handles.var.res_list {sprintf('res %d',m)}];
    end
end

set(handles.popup_sel_res,'String', handles.var.res_list);
res_act = 1;
handles.var.active_res = res_act;
set(handles.popup_sel_res,'Value', res_act);
set(handles.edit_fr,'String',num2str(handles.var.wr(res_act)/(1e9*2*pi)));
set(handles.edit_Rs,'String',num2str(handles.var.Rs(res_act)));
set(handles.edit_Q,'String',num2str(handles.var.Q(res_act)));
set(handles.slider_fr,'Value',handles.var.wr(res_act)/(1e9*2*pi));
set(handles.slider_Rs,'Value',handles.var.Rs(res_act));
set(handles.slider_Q,'Value',handles.var.Q(res_act));

handles = calcula_plota_impedancia(handles);

guidata(hObject,handles);


% --- Executes on button press in push_carrega_workspace.
function push_carrega_workspace_Callback(hObject, eventdata, handles)
% hObject    handle to push_carrega_workspace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.var.Rs = evalin('base','Rs');
handles.var.wr = evalin('base','wr');
handles.var.Q = evalin('base','Q');

handles.var.num_res = length(handles.var.Rs);

set(handles.edit_num_res,'String',num2str(handles.var.num_res));

for m=1:handles.var.num_res
    if m==1
        handles.var.res_list = {'res 1'};
    else
        handles.var.res_list = [handles.var.res_list {sprintf('res %d',m)}];
    end
end

set(handles.popup_sel_res,'String', handles.var.res_list);
res_act = 1;
handles.var.active_res = res_act;
set(handles.popup_sel_res,'Value', res_act);
set(handles.edit_fr,'String',num2str(handles.var.wr(res_act)/(1e9*2*pi)));
set(handles.edit_Rs,'String',num2str(handles.var.Rs(res_act)));
set(handles.edit_Q,'String',num2str(handles.var.Q(res_act)));
set(handles.slider_fr,'Value',handles.var.wr(res_act)/(1e9*2*pi));
set(handles.slider_Rs,'Value',handles.var.Rs(res_act));
set(handles.slider_Q,'Value',handles.var.Q(res_act));

handles = calcula_plota_impedancia(handles);

guidata(hObject,handles);


%% SALVA AJUSTE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in push_salva_em_arq.
function push_salva_em_arq_Callback(hObject, eventdata, handles)
% hObject    handle to push_salva_em_arq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
string = get(handles.popup_plano,'String');
string2 = get(handles.popup_codigo,'String');
fileName = ['Zfit_' string{handles.var.plano} '_' string2{handles.var.codigo} '.dat'];
fileName = fullfile(handles.var.path,fileName);
fp = fopen(fileName,'w+');
    fprintf(fp,'Fitting de Impedancia calculada numericamente a partir de Wake Potencial \n');
    fprintf(fp,'Plano: %s \n',string{handles.var.plano});
    fprintf(fp,'Wake Potential calculado pelo codigo: %s \n',handles.var.codigo);
    fprintf(fp,'Data: %s \n',date);
    fprintf(fp,'Rs [Ohm]    wr [Grad/s]     Q\n');
for i=1:length(handles.var.Rs)
    fprintf(fp,'%8f    %8f        %8f\n',handles.var.Rs(i),handles.var.wr(i)/1e9,handles.var.Q(i));
end
fclose(fp);

% --- Executes on button press in push_salva_workspace.
function push_salva_workspace_Callback(hObject, eventdata, handles)
% hObject    handle to push_salva_workspace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
assignin('base','Rs',handles.var.Rs);
assignin('base','wr',handles.var.wr);
assignin('base','Q',handles.var.Q);



%% Graficos e textos editados
% --- Executes during object creation, after setting all properties.
function ReZ_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ReZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



% --- Executes during object creation, after setting all properties.
function ImZ_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ImZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
