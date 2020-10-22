function varargout = EyeDetection_Gui(varargin)
% EYEDETECTION_GUI MATLAB code for EyeDetection_Gui.fig
%      EYEDETECTION_GUI, by itself, creates a new EYEDETECTION_GUI or raises the existing
%      singleton*.
%
%      H = EYEDETECTION_GUI returns the handle to a new EYEDETECTION_GUI or the handle to
%      the existing singleton*.
%
%      EYEDETECTION_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EYEDETECTION_GUI.M with the given input arguments.
%
%      EYEDETECTION_GUI('Property','Value',...) creates a new EYEDETECTION_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before EyeDetection_Gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to EyeDetection_Gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help EyeDetection_Gui

% Last Modified by GUIDE v2.5 11-Dec-2019 00:13:39

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @EyeDetection_Gui_OpeningFcn, ...
                   'gui_OutputFcn',  @EyeDetection_Gui_OutputFcn, ...
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


% --- Executes just before EyeDetection_Gui is made visible.
function EyeDetection_Gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to EyeDetection_Gui (see VARARGIN)

% Choose default command line output for EyeDetection_Gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes EyeDetection_Gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = EyeDetection_Gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in exitbutton.
function exitbutton_Callback(hObject, eventdata, handles)
close(EyeDetection_Gui);


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

sliderValue = round(get(handles.slider1,'Value'));
set(handles.edit2,'String',sliderValue);
setappdata(0,'numPart',sliderValue);

% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[200 200 200]);
end

% sliderValue = round(get(handles.slider1,'Value'));
% setappdata(0,'numPart',sliderValue);
% set(handles.edit1,'String',sliderValue)

function edit2_Callback(hObject, eventdata, handles)

editValue = round(str2num(get(handles.edit2,'String')));
set(handles.slider1,'Value',editValue);

setappdata(0,'numPart',editValue);

% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)

contents = cellstr(get(hObject,'String'));
val = contents(get(hObject,'Value'));

if(strcmp(val,'Person 2'))
    loaded_data = load('person2.csv');
    person_d = 2;
elseif(strcmp(val,'Person 5'))
    loaded_data = load('person5.csv');
    person_d = 5;
elseif(strcmp(val,'Person 6'))
    loaded_data = load('person6.csv');
    person_d = 6;
elseif(strcmp(val,'Person 7'))
    loaded_data = load('person7.csv');
    person_d = 7;
end

setappdata(0,'person_d',person_d);
setappdata(0,'loaded_data',loaded_data);

% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit1_Callback(hObject, eventdata, handles)

editValue = round(str2num(get(handles.edit1,'String')));
setappdata(0,'process_noise',editValue);

% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)

editValue = round(str2num(get(handles.edit3,'String')));
setappdata(0,'measurement_noise',editValue);


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in start.
function start_Callback(hObject, eventdata, handles)
% hObject    handle to start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

person_d = getappdata(0,'person_d');
measurement_noise = getappdata(0,'measurement_noise');
process_noise = getappdata(0,'process_noise');
numPart = getappdata(0,'numPart');
person = getappdata(0,'loaded_data');
accuracy = getappdata(0,'accuracy');

video = VideoWriter(['person' num2str(person_d) '.avi']); %create the video object

video.FrameRate = 5;
open(video); %open the file for writing

img = imread(['00' num2str(person_d) '_01.png']);
img = rgb2gray(img);

% Read the Data for eye's position
LeftEyeCenter = person(:,3:4);
RightEyeCenter = person(:,9:10);

% Variables for marking the ground truth
r = 3;
th = 0:pi/50:2*pi;
k = 5; % Number of Observations per image

% Particle Filter parameters
params.Sigma_Q = diag([measurement_noise measurement_noise]); % measurement noise covariance matrix
params.M = numPart;                    % Number of particles to use
params.state_space_bound = [size(img,1);size(img,2)];
params.Sigma_R = diag([process_noise process_noise]);   % process noise covariance matrix
params.thresh_avg_likelihood = 0.0001;
params.v_0 = 2*pi*200/k;
params.theta_0 = 0;

% Variable Initialization
S.X = [rand(1, params.M)*params.state_space_bound(1); % sampling uniformly from the state space
     rand(1, params.M)*params.state_space_bound(2)];
mean_S = zeros(2,k);
S.W = 1/params.M * ones(1,params.M); % initialize equal weights

S_R.X = [rand(1, params.M)*params.state_space_bound(1); % sampling uniformly from the state space
     rand(1, params.M)*params.state_space_bound(2)];
 
mean_R_S = zeros(2,k);
S_R.W = 1/params.M * ones(1,params.M); % initialize equal weights

X_L = [];
Z_L = [];

X_R = [];
Z_R = [];

mean_left = [];
mean_right = [];

axes(handles.axes1);
for d = 1:3
    for i = 1:size(person,1)
        
        % Read the 24 images
        if i < 10
            c = ['00' num2str(person_d) '_0' num2str(i) '.png'];
        elseif i < 13
            c = ['00' num2str(person_d) '_' num2str(i) '.png'];
        elseif i < 16
            c = ['00' num2str(person_d) '_' num2str(24-i+1) '.png']; 
        else
            c = ['00' num2str(person_d) '_0' num2str(24-i+1) '.png'];
        end      

        img = imread(c); 

        measurements_Left = (accuracy * (rand(2,k) - 0.5)) + LeftEyeCenter(i,:)';
        measurements_Right = (accuracy * (rand(2,k) - 0.5)) + RightEyeCenter(i,:)';

        X_L = [X_L repmat(LeftEyeCenter(i,:)',1,k)];
        Z_L = [Z_L measurements_Left];
        
        X_R = [X_R repmat(RightEyeCenter(i,:)',1,k)];
        Z_R = [Z_R measurements_Right];
       
        % Create Circles around the eyes for each image
        circle_left = ([r * cos(th) + LeftEyeCenter(i,1);r * sin(th) + LeftEyeCenter(i,2)]);
        circle_right = ([r * cos(th) + RightEyeCenter(i,1);r * sin(th) + RightEyeCenter(i,2)]);

        % To create red boxes on the person's eye
        imnew = zeros(size(img,1),size(img,2),3);
        imnew(round(circle_left(2,:)),round(circle_left(1,:)),2) = 255;
        imnew(round(circle_right(2,:)),round(circle_right(1,:)),2) = 255;

        % Combine the two images
        image = img + uint8(imnew);
        imshow(image);
        
        %%%%%%% Particle Filter Algorithm Implementation
        for j = 1 : k               
            %%%%% Left Eye
            S_bar = pf_predict(S,params);
            z_L = measurements_Left;
            S_bar = pf_weight(S_bar,z_L,params);
            S = systematic_resample(S_bar);

            %%%%% Right Eye
            S_bar_R = pf_predict(S_R,params);
            z_R = measurements_Right;
            S_bar_R = pf_weight(S_bar_R,z_R,params);
            S_R = systematic_resample(S_bar_R);

            mean_S(:,j) = mean(S.X(1:2,:),2); % compute the center of the distribution
            mean_R_S(:,j) = mean(S_R.X(1:2,:),2); % compute the center of the distribution

            hold on;
            plot(z_L(1,j),z_L(2,j),'rx','MarkerSize',10);      
            plot(z_R(1,j),z_R(2,j),'rx','MarkerSize',10);

            plot(S.X(1,:),S.X(2,:),'b.');
            plot(S_R.X(1,:),S_R.X(2,:),'b.');
            hold off;

            drawnow;
        end
        %%%%%%%

        writeVideo(video,image); % Add image to the video
        mean_left = [mean_left mean_S];
        mean_right = [mean_right mean_R_S];
               
    end
end
close(video); % Close the video and save it

err_z = sqrt(sum((Z_L - X_L).^2,1));
err_xhat = sqrt(sum((X_L - mean_left).^2,1));
merr_z = mean(err_z);
merr_xhat = mean(err_xhat);
stderr_z = std(err_z);
stderr_xhat = std(err_xhat);
format compact;

figure('Name','Left Eye');
clf;
plot(err_z,'r-');
hold on;
plot(err_xhat,'b-');
title(sprintf('absolute error analysis: measurements: %0.1f \\pm %0.1f',merr_z,stderr_z));

err_z = sqrt(sum((Z_R - X_R).^2,1));
err_xhat = sqrt(sum((X_R - mean_right).^2,1));
merr_z = mean(err_z);
merr_xhat = mean(err_xhat);
stderr_z = std(err_z);
stderr_xhat = std(err_xhat);
format compact;

figure('Name','Right Eye');
clf;
plot(err_z,'r-');
hold on;
plot(err_xhat,'b-');
title(sprintf('absolute error analysis: measurements: %0.1f \\pm %0.1f',merr_z,stderr_z));



function edit4_Callback(hObject, eventdata, handles)

editValue = round(str2num(get(handles.edit4,'String')));
setappdata(0,'accuracy',editValue);

% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
