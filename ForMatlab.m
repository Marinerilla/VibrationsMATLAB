function varargout = ForMatlab(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%                                            %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%  The following script is divised by:       %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%                                            %%%%%%%%%%%                           
%%%%%%%%%%%%%%%%%%%%    Name :                                  %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%    Module :                                %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%                                            %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%                                            %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last Modified by GUIDE v2.5 18-Apr-2013 18:38:10

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ForMatlab_OpeningFcn, ...
                   'gui_OutputFcn',  @ForMatlab_OutputFcn, ...
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


% --- Executes just before ForMatlab is made visible.
function ForMatlab_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ForMatlab (see VARARGIN)

% Choose default command line output for ForMatlab
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

%%initialize_gui(hObject, handles, false);

% UIWAIT makes ForMatlab wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ForMatlab_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function length_Callback(hObject, eventdata, handles)

%        get(hObject,'String') returns contents of length as text
%        str2double(get(hObject,'String')) returns contents of length as a double
length = str2double(get(hObject, 'String'));
if isnan(length)
    set(hObject, 'String', '');
    errordlg('Input must be a number','Error');
end

% Save the new length value
handles.metricdata.length = length;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function length_CreateFcn(hObject, eventdata, handles)


% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function radius_Callback(hObject, eventdata, handles)


%        get(hObject,'String') returns contents of radius as text
%        str2double(get(hObject,'String')) returns contents of radius as a double
radius = str2double(get(hObject, 'String'));
if isnan(radius)
    set(hObject, 'String', '');
    errordlg('Input must be a number','Error');
end

% Save the new radius value
handles.metricdata.radius = radius;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function radius_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function youngm_Callback(hObject, eventdata, handles)

%        get(hObject,'String') returns contents of youngm as text
%        str2double(get(hObject,'String')) returns contents of youngm as a double
youngm = str2double(get(hObject, 'String'));
if isnan(youngm)
    set(hObject, 'String', '');
    errordlg('Input must be a number','Error');
end

% Save the new youngm value
handles.metricdata.youngm = youngm;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function youngm_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function density_Callback(hObject, eventdata, handles)

%        get(hObject,'String') returns contents of density as text
%        str2double(get(hObject,'String')) returns contents of density as a double
density = str2double(get(hObject, 'String'));
if isnan(density)
    set(hObject, 'String', '');
    errordlg('Input must be a number','Error');
end

% Save the new density value
handles.metricdata.density = density;
guidata(hObject,handles)

function density_CreateFcn(hObject, eventdata, handles)


% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function dens_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in bc_popupmenu.
function bc_popupmenu_Callback(hObject, eventdata, handles)



% --- Executes during object creation, after setting all properties.
function bc_popupmenu_CreateFcn(hObject, eventdata, handles)


% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in evaluate_pushbutton.
function evaluate_pushbutton_Callback(hObject, eventdata, handles)

% makes L,R,E,Ro handle the metricdata entered by the user in the input
L = handles.metricdata.length;
R = handles.metricdata.radius;
E = handles.metricdata.youngm;
Ro = handles.metricdata.density;

%Calculates moment of Inertia and Area based on radius
Ix=(1/4)*pi*R^4;
A=pi*R^2;

% Initial guess values for fzero to obtain the first five natural frequency and
% modeshapes
aO = [4.7 7.2 10.8 14 17];
zO = [1.8 4.7 7.2 10.8 14];
uO = [1.5 3.7 4.7 7 7.7];
pO = [3 6.2 9.3 12.5 15.6];

%Creates a loop for i, from first value to fifth value in steps of 1
for i=1:1:5
    
% Obtains the spatial frecuencies using fzero with the initial guess values
an(i)=fzero('cos(x)*cosh(x)-1',aO(i));
zn(i)=fzero('cos(x)*cosh(x)+1',zO(i));
un(i)=fzero('tan(x)-tanh(x)',uO(i));
pn(i)=fzero('sin(x)',pO(i));

% coeficients Qn(L)/Sn(L) and Sn(L)/Pn(L) theres no need for coeficients
% for Pinned-Pinned boundary condition
F(i)=(cosh(an(i))-cos(an(i)))./(sinh(an(i))-sin(an(i)));    % Qn(L)/Sn(L) for Free-Free & Fixed-Fixed
K(i)=(sinh(zn(i))-sin(zn(i)))./(cosh(zn(i))-cos(zn(i)));    % Sn(L)/Pn(L) for Fixed-Free
Q(i)=(cosh(un(i))-cos(un(i)))./(sinh(un(i))-sin(un(i)));    % Qn(L)/Sn(L) for Fixed-Pinned

% Divides the spatial frecuencies by L to get the Bn value for each
% boundary condition
ax(i)=(an(i)/L)';
zx(i)=(zn(i)/L)';
ux(i)=(un(i)/L)';
px(i)=(pn(i)/L)';

% natural frecuency Wn = 2*pi*fn 
% fn = Wn/(2*pi)=(Bn^2*(SQRT((E*I)/(Ro*A)))/2*pi
% fn = natural frecuency
fna(i) =(ax(i)^2*sqrt((E*Ix)/(Ro*A)))/(2*pi); % Natural Frequency for Free-Free & Fixed-Fixed Boundary Condition
fnz(i)=(zx(i)^2*sqrt((E*Ix)/(Ro*A)))/(2*pi);  % Natural Frequency for Fixed-Free Boundary Condition
fnu(i)=(ux(i)^2*sqrt((E*Ix)/(Ro*A)))/(2*pi);  % Natural Frequency for Fixed-Pinned Boundary Condition
fnp(i)=(px(i)^2*sqrt((E*Ix)/(Ro*A)))/(2*pi);  % Natural Frequency for Pinned-Pinned Boundary Condition

end

% Creates a vector of 180 linearly spaced numbers from 0 to L
x=linspace(0, L, 180);
xl=x./L; %Makes certain that the boundary conditions of the beam corresponds to that of the modeshapes

Tc0='(cosh(ax(i).*x(j))+cos(ax(i).*x(j)))-F(i).*(sinh(ax(i).*x(j))+sin(ax(i).*x(j)))';  % Free-Free 
Tc1='(cosh(zx(i).*x(j))-cos(zx(i).*x(j)))-K(i).*(sinh(zx(i).*x(j))-sin(zx(i).*x(j)))';  % Fixed-Fixed
Tc2='(cosh(ux(i).*x(j))-cos(ux(i).*x(j)))-Q(i).*(sinh(ux(i).*x(j))-sin(ux(i).*x(j)))';  % Fixed-Free
Tc3='(cosh(ax(i).*x(j))-cos(ax(i).*x(j)))-F(i).*(sinh(ax(i).*x(j))-sin(ax(i).*x(j)))';  % Fixed-Pinned
Tc4='(sin(px(i)*x(j)))'; % Pinned-Pinned

% Generates an array of zero's for all the coloums x. This is for the
% boundary conditions and spartial frequencies
Xnx0=zeros(5,length(x)); % Free-Free 
Xnx1=zeros(5,length(x)); % Fixed-Fixed
Xnx2=zeros(5,length(x)); % Fixed-Free
Xnx3=zeros(5,length(x)); % Fixed-Pinned
Xnx4=zeros(5,length(x)); % Pinned-Pinned

% loop to obtain all the values of each boundary condition mode shape,
% defining the matrices, the components of each matrix row must be the
% solution for each "x" from the linspace of the beam of the equation of
% each mode shape of each boundary condition
for i=1:1:5
    for j=1:length(x)
        Xnx0(i,j)=eval(Tc0);
        Xnx1(i,j)=eval(Tc1);
        Xnx2(i,j)=eval(Tc2);
        Xnx3(i,j)=eval(Tc3);
        Xnx4(i,j)=eval(Tc4);% Pinned-Pinned
    end

end


axes(handles.axes1);
cla;

% Obtain the boundary conditions
popup_sel_index = get(handles.bc_popupmenu, 'Value');

% Generates plot based on the case selected, i.e. based on the boundary
% conditions selected
switch popup_sel_index
    case 1   % Free-Free
         handles.text.n1 = fna(i);
    handles.text.n2 = fna(i);
    handles.text.n3 = fna(i);
    handles.text.n4 = fna(i);
    handles.text.n5 = fna(i);
    set(handles.n1, 'String', fna(1,1))
    set(handles.n2, 'String', fna(1,2))
    set(handles.n3, 'String', fna(1,3))
    set(handles.n4, 'String', fna(1,4))
    set(handles.n5, 'String', fna(1,5))
         plot(xl,Xnx0(1,:), 'b-');hold on;
    plot(xl,Xnx0(2,:), 'r-');
    plot(xl,Xnx0(3,:), 'g-');
    plot(xl,Xnx0(4,:), 'm-');
    plot(xl,Xnx0(5,:), 'k-');
    title('Mode shape of the Free-free');
    legend('Mode #1','Mode #2','Mode #3','Mode #4','Mode #5', 0); xlabel('x/L'); ylabel('Mode shape'); grid;
    
    
    case 2 % Fixed-Free
         handles.text.n1 = fnz(i);
    handles.text.n2 = fnz(i);
    handles.text.n3 = fnz(i);
    handles.text.n4 = fnz(i);
    handles.text.n5 = fnz(i);
    set(handles.n1, 'String', fnz(1,1));
    set(handles.n2, 'String', fnz(1,2));
    set(handles.n3, 'String', fnz(1,3));
    set(handles.n4, 'String', fnz(1,4));
    set(handles.n5, 'String', fnz(1,5));
         plot(xl,Xnx1(1,:), 'b-');hold on;
    plot(xl,Xnx1(2,:), 'r-');
    plot(xl,Xnx1(3,:), 'g-');
    plot(xl,Xnx1(4,:), 'm-');
    plot(xl,Xnx1(5,:), 'k-');
    title('Mode shape of the Fixed-free');
    legend('Mode #1','Mode #2','Mode #3','Mode #4','Mode #5', 0); xlabel('x/L'); ylabel('Mode shape'); grid;  
    
    case 3 % Fixed-Pinned
         handles.text.n1 = fna(i);
    handles.text.n2 = fnu(i);
    handles.text.n3 = fnu(i);
    handles.text.n4 = fnu(i);
    handles.text.n5 = fnu(i);
    set(handles.n1, 'String', fnu(1,1))
    set(handles.n2, 'String', fnu(1,2))
    set(handles.n3, 'String', fnu(1,3))
    set(handles.n4, 'String', fnu(1,4))
    set(handles.n5, 'String', fnu(1,5))
         plot(xl,Xnx2(1,:), 'b-');hold on;
    plot(xl,Xnx2(2,:), 'r-');
    plot(xl,Xnx2(3,:), 'g-');
    plot(xl,Xnx2(4,:), 'm-');
    plot(xl,Xnx2(5,:), 'k-');
    title('Mode shape of the Fixed-pinned');
    legend('Mode #1','Mode #2','Mode #3','Mode #4','Mode #5', 0); xlabel('x/L'); ylabel('Mode shape'); grid;
    
    case 4 % Fixed-Fixed
        handles.text.n1 = fna(i);
    handles.text.n2 = fna(i);
    handles.text.n3 = fna(i);
    handles.text.n4 = fna(i);
    handles.text.n5 = fna(i);
    set(handles.n1, 'String', fna(1,1))
    set(handles.n2, 'String', fna(1,2))
    set(handles.n3, 'String', fna(1,3))
    set(handles.n4, 'String', fna(1,4))
    set(handles.n5, 'String', fna(1,5))
         plot(xl,Xnx3(1,:), 'b-');hold on;
    plot(xl,Xnx3(2,:), 'r-');
    plot(xl,Xnx3(3,:), 'g-');
    plot(xl,Xnx3(4,:), 'm-');
    plot(xl,Xnx3(5,:), 'k-');
    title('Mode shape of the Fixed-fixed');
    legend('Mode #1','Mode #2','Mode #3','Mode #4','Mode #5', 0); xlabel('x/L'); ylabel('Mode shape'); grid;
    
    case 5 % Pinned-Pinned
        handles.text.n1 = fna(i);
    handles.text.n2 = fnp(i);
    handles.text.n3 = fnp(i);
    handles.text.n4 = fnp(i);
    handles.text.n5 = fnp(i);
    set(handles.n1, 'String', fnp(1,1))
    set(handles.n2, 'String', fnp(1,2))
    set(handles.n3, 'String', fnp(1,3))
    set(handles.n4, 'String', fnp(1,4))
    set(handles.n5, 'String', fnp(1,5))
         plot(xl,Xnx4(1,:), 'b-');hold on;
    plot(xl,Xnx4(2,:), 'r-');
    plot(xl,Xnx4(3,:), 'g-');
    plot(xl,Xnx4(4,:), 'm-');
    plot(xl,Xnx4(5,:), 'k-');
    title('Mode shape of the Pinned-pinned');
    legend('Mode #1','Mode #2','Mode #3','Mode #4','Mode #5', 0); xlabel('x/L'); ylabel('Mode shape'); grid;
    
end


% --- Executes on button press in animate_pushbutton.
function animate_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to animate_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in reset_pushbutton.
function reset_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to reset_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
initialize_gui(gcbf, handles, true);



% --------------------------------------------------------------------
function initialize_gui(fig_handle, handles, isreset)
% If the metricdata field is present and the reset flag is false, it means
% we are we are just re-initializing a GUI by calling it from the cmd line
% while it is up. So, bail out as we dont want to reset the data.
if isfield(handles, 'metricdata') && ~isreset
    return;
end

% The following codes tell matlab that if the reset flag is true, then
% handle and set all the input fields to an empty value and the text to 0
handles.metricdata.length  ='';
handles.metricdata.radius  ='';
handles.metricdata.youngm  ='';
handles.metricdata.density ='';

handles.text.n1 ='';
handles.text.n2 ='';
handles.text.n3 ='';
handles.text.n4 ='';
handles.text.n5 ='';

set(handles.density, 'String', handles.metricdata.density);
set(handles.length,  'String', handles.metricdata.length);
set(handles.radius,  'String', handles.metricdata.radius);
set(handles.youngm, 'String', '');

set(handles.n1, 'String', '0');
set(handles.n2, 'String', '0');
set(handles.n3, 'String', '0');
set(handles.n4, 'String', '0');
set(handles.n5, 'String', '0');


% Update handles structure
guidata(handles.figure1, handles);


function n1_CreateFcn(hObject, eventdata, handles)




function tn1_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function tn1_CreateFcn(hObject, eventdata, handles)


% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
