%% Author: Abdulellah Abualshour
%  King Abdullah University of Science and Technology

function varargout = rec_gui(varargin)
% REC_GUI MATLAB code for rec_gui.fig
%      REC_GUI, by itself, creates a new REC_GUI or raises the existing
%      singleton*.
%
%      H = REC_GUI returns the handle to a new REC_GUI or the handle to
%      the existing singleton*.
%
%      REC_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in REC_GUI.M with the given input arguments.
%
%      REC_GUI('Property','Value',...) creates a new REC_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before rec_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to rec_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help rec_gui

% Last Modified by GUIDE v2.5 10-Jul-2017 10:48:14

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @rec_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @rec_gui_OutputFcn, ...
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


% --- Executes just before rec_gui is made visible.
function rec_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to rec_gui (see VARARGIN)

% Choose default command line output for rec_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes rec_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = rec_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename pathname] = uigetfile({'*.jpg';'*.bmp'},'File Selector');
 handles.myImage = strcat(pathname, filename);
 display(handles.myImage);
 axes(handles.axes1);
 imshow(handles.myImage);
 handles.A = imread(handles.myImage);
 handles.Backup = handles.A;
 % save the updated handles object
 guidata(hObject,handles);


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
xoffset = str2double(get(handles.edit3,'String'));
yoffset = str2double(get(handles.edit4,'String'));
display(xoffset);
display(yoffset);

B = imtranslate(handles.A,[xoffset,yoffset],'OutputView','full');
handles.A = B;
imshow(handles.A);
% save the updated handles object
guidata(hObject,handles);




% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
angle = str2double(get(handles.edit2,'String'));
display(angle);
B = imrotate(handles.A, angle);
handles.A = B;
imshow(handles.A);
% save the updated handles object
guidata(hObject,handles);




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


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
I = rgb2gray(handles.A);

L = lsd(double(I));
Lall = L;
angle_threshold=20;

%% VERTICAL PERSPECTIVE
%segment endpoints
X1=L(:,1); X2=L(:,3); Y1=L(:,2); Y2=L(:,4);
XY=[X2-X1 Y2-Y1];

%segment orientations
lineangle =(atan(XY(:,2)./XY(:,1)))/pi*180;
lineangle(lineangle<0)=lineangle(lineangle<0)+180;

L(lineangle<angle_threshold | lineangle>(180-angle_threshold),:)=[];

midpoints = (L(:,1:2)+L(:,3:4))/2;

%segment endpoints
X1=L(:,1); X2=L(:,3); Y1=L(:,2); Y2=L(:,4);
XY=[X2-X1 Y2-Y1];

%segment orientations
lineangle =(atan(XY(:,2)./XY(:,1)))/pi*180;
lineangle(lineangle<0)=lineangle(lineangle<0)+180;

data = [midpoints(:,1)  lineangle]';

iterations = 20; % k
min_points = 2; % n
threshold = 4; % t
num_inliers_required = 140; % d

% initialization
num_inliers_best = 0;
x = data(1,:);
y = data(2,:);

% main loop
for i=1:iterations
    
    % randomly select points from data (subset)
    [~,idx] = datasample(x,min_points);
    % get coffecients
    coef = polyfit(x(1,idx),y(1,idx),1); % uses least squares
    % get predictions of y from coefficients using polyval
    y_pred = polyval(coef,x);
    % find square residual
    square_resid = (y - y_pred).^2;
    % compare all points (using square residual) to threshold
    inliers = [x(square_resid < threshold); y(square_resid < threshold)];
    inliers_x1 = X1(square_resid < threshold);
    inliers_x2 = X2(square_resid < threshold);
    inliers_y1 = Y1(square_resid < threshold);
    inliers_y2 = Y2(square_resid < threshold);
    % skip iteration if size of current inliers < size of best
    if length(inliers) < num_inliers_best
        continue;
    end
    % current random sample is best sample if reached here
    best_inliers = inliers;
    
    num_inliers_best = length(best_inliers);
    best_model = coef;
    % break out of the loop if we reached the required # of inliers
    if num_inliers_best >= num_inliers_required
        break
    end
    
end

% optimization:
a = [inliers_x1 inliers_y1];

b = [inliers_x2 inliers_y2];

imageSize = size(I);
imageWidth = imageSize(1);
imageLength = imageSize(2);

%function
%fun = @(a,b)symsum(abs(((h1'*a(i))/(h3'*a(i))) - ((h1'*a(i))/(h3'*a(i)))), i, 1, k);

fun = @(d) rect_func_v(d,imageLength,imageWidth,a',b');

[d,f] = fminsearch(fun,[0 0]);
 
[~,Hv] = rect_func_v(d,imageLength,imageWidth,[],[]);

boundaries=Hv*[1 imageWidth imageWidth 1;1 1 imageLength imageLength;1 1 1 1];
boundaries=boundaries./repmat(boundaries(3,:),[3 1]);
x=boundaries(1,1);

width_v=boundaries(1,2)+x;
length_v=boundaries(2,3);


%transform original lines 
L = Lall;
tform = maketform('projective',Hv'); % 2d spatial transformation required
% transforming image to affine using H:
newImage = imtransform(handles.A, tform);

A=L(:,1:2)';B=L(:,3:4)';k=size(A,2);
A_=Hv*[A;ones(1,k)];A_=A_./repmat(A_(3,:),[3 1]);
B_=Hv*[B;ones(1,k)];B_=B_./repmat(B_(3,:),[3 1]);
Lv=[A_(1,:)-x;A_(2,:);B_(1,:)-x;B_(2,:)]';


%% HORIZONTAL PERSPECTIVE

angle_threshold=20;

%segment endpoints
X1=Lv(:,1); X2=Lv(:,3); Y1=Lv(:,2); Y2=Lv(:,4);
XY=[X2-X1 Y2-Y1];

%segment orientations
lineangle =(atan(XY(:,2)./XY(:,1)))/pi*180;

Lv(abs(lineangle)>(90-angle_threshold),:)=[];

midpoints = (Lv(:,1:2)+Lv(:,3:4))/2;

%segment endpoints
X1=Lv(:,1); X2=Lv(:,3); Y1=Lv(:,2); Y2=Lv(:,4);
XY=[X2-X1 Y2-Y1];

%segment orientations
lineangle =(atan(XY(:,2)./XY(:,1)))/pi*180;
%lineangle(lineangle<0)=lineangle(lineangle<0)+180;

data = [midpoints(:,2)  lineangle]';

iterations = 20; % k
min_points = 2; % n
threshold = 7; % t
num_inliers_required_h = 140; % d

% initialization
num_inliers_best_h = 0;
x = data(1,:);
y = data(2,:);

% main loop
for i=1:iterations
    
    % randomly select points from data (subset)
    [~,idx] = datasample(x,min_points);
    % get coffecients
    coef = polyfit(x(1,idx),y(1,idx),1); % uses least squares
    % get predictions of y from coefficients using polyval
    y_pred = polyval(coef,x);
    % find square residual
    square_resid = (y - y_pred).^2;
    % compare all points (using square residual) to threshold
    inliers_h = [x(square_resid < threshold); y(square_resid < threshold)];
    inliers_x1_h = X1(square_resid < threshold);
    inliers_x2_h = X2(square_resid < threshold);
    inliers_y1_h = Y1(square_resid < threshold);
    inliers_y2_h = Y2(square_resid < threshold);
    % skip iteration if size of current inliers < size of best
    if length(inliers_h) < num_inliers_best_h
        continue;
    end
    % current random sample is best sample if reached here
    best_inliers_h = inliers_h;
    
    num_inliers_best_h = length(best_inliers_h);
    best_model = coef;
    % break out of the loop if we reached the required # of inliers
    if num_inliers_best_h >= num_inliers_required_h
        break
    end
    
end

% optimization:
a = [inliers_x1_h inliers_y1_h];

b = [inliers_x2_h inliers_y2_h];

fun = @(d) rect_func_h(d,length_v,width_v,a',b');

options = optimset;
options = optimset(options,'Display', 'off');
options = optimset(options,'Algorithm', 'active-set');

[d,f] = fmincon(fun,[0 0],[],[],[],[],[-4*imageLength -4*imageLength],[4*imageLength 4*imageLength],[],options);

[~,Hh] = rect_func_h(d,imageLength,imageWidth,[],[]);


%transform
tform = maketform('projective',Hh'); % 2d spatial transformation required
% transforming image to affine using H:
finalImage = imtransform(newImage, tform);

imshow(finalImage);




% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



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


% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[x,y] = snake_read;
handles.p1 = [x(1) y(1) 1];
handles.p2 = [x(2) y(2) 1];
handles.p3 = [x(3) y(3) 1];
handles.p4 = [x(4) y(4) 1];
handles.x = x;
handles.y = y;
guidata(hObject,handles);




% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.A = handles.Backup;
imshow(handles.A);
% save the updated handles object
guidata(hObject,handles);


% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% locate the pair of parallel lines:
l1 = cross(handles.p1,handles.p2)
l2 = cross(handles.p3,handles.p4);
l3 = cross(handles.p2, handles.p3);
l4 = cross(handles.p1, handles.p4);
% find the  line at infinity:
infA = cross(l1, l2);
infA = infA/infA(1,3); % first pair
infB = cross(l3, l4);
infB = infB/ infB(1,3); % second pair
l = cross(infA, infB) % line at infinity
% transformation matrix H, put homogeneous coordinates of line at infinity:
l(1, 1);
l(1, 2);
l(1, 3);
H = [1 0 0; 0 1 0; l(1, 1)/l(1,3) l(1, 2)/l(1,3) 1];
tform = maketform('projective',H'); % 2d spatial transformation required
% transforming image to affine using H:
affineImage = imtransform(handles.A, tform);
handles.p1'
p1x = H* handles.p1'
p1x = [p1x(1)/p1x(3); p1x(2)/p1x(3); 1]
%p1x = p1x';
p2x = H* handles.p2'
p2x = [p2x(1)/p2x(3); p2x(2)/p2x(3); 1]
%p2x = p2x';
p3x = H* handles.p3'
p3x = [p3x(1)/p3x(3); p3x(2)/p3x(3); 1]
%p3x = p3x';
p4x = H* handles.p4'
p4x = [p4x(1)/p4x(3); p4x(2)/p4x(3); 1]
%p4x = p4x';

% begin metric transformation
% find l_a and l_b, define theta

la_1 = cross(p1x,p2x);
la_1 = [la_1(1)/la_1(3) la_1(2)/la_1(3) 1]

lb_1 = cross(p4x,p1x);
lb_1 = [lb_1(1)/lb_1(3) lb_1(2)/lb_1(3) 1]

la_2 = cross(p2x,p3x);
la_2 = [la_2(1)/la_2(3) la_2(2)/la_2(3) 1]

lb_2 = cross(p3x,p4x);
lb_2 = [lb_2(1)/lb_2(3) lb_2(2)/lb_2(3) 1]

% compute variables a and b:
a1 = -la_1(2)/la_1(1)
a2 = -la_2(2)/la_2(1)
b1 = -lb_1(2)/lb_1(1)
b2 = -lb_2(2)/lb_2(1)
% find radius and centers of circle
r1 = abs(a1 - b1)/2
r2 = abs(a2 - b2)/2
c1_alpha = (a1 + b2)/2
c2_alpha = (a2 + b1)/2
c1_beta = 0;
c2_beta = 0;
% intersection of the two circles:
[pt1, pt2] = circcirc(c1_alpha,c1_beta,r1,c2_alpha,c2_beta,r2)
% pt1 is the point we need
alpha = pt1(1);
beta = pt1(2);
% define transformation matrix:
H = [(1/beta) -(alpha/beta) 0; 0 1 0; 0 0 1];
tform = maketform('projective',H');
metricImage = imtransform(affineImage, tform);

imshow(metricImage);