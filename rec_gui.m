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
% segment endpoints
X1=L(:,1); X2=L(:,3); Y1=L(:,2); Y2=L(:,4);
XY=[X2-X1 Y2-Y1];

% segment orientations
lineangle =(atan(XY(:,2)./XY(:,1)))/pi*180;
lineangle(lineangle<0)=lineangle(lineangle<0)+180;
L(lineangle<angle_threshold | lineangle>(180-angle_threshold),:)=[];
midpoints = (L(:,1:2)+L(:,3:4))/2;

% segment endpoints
X1=L(:,1); X2=L(:,3); Y1=L(:,2); Y2=L(:,4);
XY=[X2-X1 Y2-Y1];

% segment orientations
lineangle =(atan(XY(:,2)./XY(:,1)))/pi*180;
lineangle(lineangle<0)=lineangle(lineangle<0)+180;
data = [midpoints(:,1)  lineangle]';

% RANSAC
iterations = 20; % k
min_points = 2; % n
threshold = 4; % t
num_inliers_required = 140; % d
num_inliers_best = 0;
x = data(1,:);
y = data(2,:);

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

% function
fun = @(d) rect_func_v(d,imageLength,imageWidth,a',b');
[d,f] = fminsearch(fun,[0 0]);
[~,Hv] = rect_func_v(d,imageLength,imageWidth,[],[]);
boundaries=Hv*[1 imageWidth imageWidth 1;1 1 imageLength imageLength;1 1 1 1];
boundaries=boundaries./repmat(boundaries(3,:),[3 1]);
x=boundaries(1,1);
width_v=boundaries(1,2)+x;
length_v=boundaries(2,3);

% transform original lines 
L = Lall;
tform = maketform('projective',Hv'); % 2d spatial transformation required

% transforming image to affine using H:
newImage = imtransform(handles.A, tform);

% prepare Lv for horizontal perspective
A=L(:,1:2)';B=L(:,3:4)';k=size(A,2);
A_=Hv*[A;ones(1,k)];A_=A_./repmat(A_(3,:),[3 1]);
B_=Hv*[B;ones(1,k)];B_=B_./repmat(B_(3,:),[3 1]);
Lv=[A_(1,:)-x;A_(2,:);B_(1,:)-x;B_(2,:)]';


%% HORIZONTAL PERSPECTIVE

angle_threshold=20;

% segment endpoints
X1=Lv(:,1); X2=Lv(:,3); Y1=Lv(:,2); Y2=Lv(:,4);
XY=[X2-X1 Y2-Y1];

% segment orientations
lineangle =(atan(XY(:,2)./XY(:,1)))/pi*180;
Lv(abs(lineangle)>(90-angle_threshold),:)=[];
midpoints = (Lv(:,1:2)+Lv(:,3:4))/2;

% segment endpoints
X1=Lv(:,1); X2=Lv(:,3); Y1=Lv(:,2); Y2=Lv(:,4);
XY=[X2-X1 Y2-Y1];

% segment orientations
lineangle =(atan(XY(:,2)./XY(:,1)))/pi*180;


% RANSAC
data = [midpoints(:,2)  lineangle]';
iterations = 20; % k
min_points = 2; % n
threshold = 7; % t
num_inliers_required_h = 140; % d
num_inliers_best_h = 0;
x = data(1,:);
y = data(2,:);

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

% transform
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
im = handles.A;
p1 = handles.p1;
p2 = handles.p2;
p3 = handles.p3;
p4 = handles.p4;

l1 = cross(p1,p2)
l2 = cross(p3,p4);
l3 = cross(p2, p3);
l4 = cross(p1, p4);

% find the  line at infinity:
a = cross(l1, l2);
a = a/a(1,3); % first pair
b = cross(l3, l4);
b = b/ b(1,3); % second pair

l_inf = cross(a, b);
l_inf=l_inf/l_inf(3);

% Recover the projective transformation that maps the line at infinity to its canonical position
P = [1,0,0;0,1,0;l_inf(1),l_inf(2),1];

m1 = l4;
m2 = l3;
v1=cross(l1,l2);v1=v1/v1(3);
v2=cross(m1,m2);v2=v2/v2(3);

% Define constraints
% contraint: a known right angle between two orthogonal lines
p1x = P* handles.p1'
p1x = [p1x(1)/p1x(3); p1x(2)/p1x(3); 1]
p2x = P* handles.p2'
p2x = [p2x(1)/p2x(3); p2x(2)/p2x(3); 1]
p3x = P* handles.p3'
p3x = [p3x(1)/p3x(3); p3x(2)/p3x(3); 1]
p4x = P* handles.p4'
p4x = [p4x(1)/p4x(3); p4x(2)/p4x(3); 1]

la = cross(p1x,p2x);
la = [la(1)/la(3) la(2)/la(3) 1]
lb = cross(p2x,p3x);
lb = [lb(1)/lb(3) lb(2)/lb(3) 1]
a1 = -la(2)/la(1)
b1 = -lb(2)/lb(1)
rd1 = abs(a1 - b1)/2
c1_alpha = (a1 + b1)/2
c1_beta = 0;

% contraint: angle between diagonals given known length ratio
% calculate length ratio
inter1 = cross(l1,m1); p1 = inter1/inter1(3);
inter2 = cross(l2,m1); p2 = inter2/inter2(3);
inter3 = cross(l1,m2); p3 = inter3/inter3(3);
inter4 = cross(l2,m2); p4 = inter4/inter4(3);

% estimate aspect ratio through rectangle vertices
k2=dot((cross(p1,p4)),p3)/dot((cross(p2,p4)),p3);
k3=dot((cross(p1,p4)),p2)/dot((cross(p3,p4)),p2);
n2(:,1)=k2*p2-p1;
n3(:,1)=k3*p3-p1;
im_width=size(im,2);
im_height=size(im,1);
u0=im_width/2;
v0=im_height/2;
f=(-v1(1)*u0-v1(2)*v0-1)/(v1(1)*v1(2)+v2(1)*v2(2));
A=[f 0 u0;0 f v0;0 0 1];
s=(n2'*inv(A)'*(A\n2))/...
    (n3'*inv(A)'*(A\n3));
s=sqrt(s);

% apply ratio to metric constraints
ld1=cross(p1,p4);ld1=ld1/ld1(3);ld1A=ld1/P;
ld2=cross(p2,p3);ld2=ld2/ld2(3);ld2A=ld2/P;
theta=2*asin(1/sqrt(s^2+1));
a2=-ld1A(2)/ld1A(1);b2=-ld2A(2)/ld2A(1);
c2_alpha=(a2+b2)/2;
c2_beta=((a2-b2)/2)*cot(theta);
rd2=abs((a2-b2)/(2*sin(theta)));

% Constraint combination
if rd1>0 && rd2>0 &&rd1~=Inf && rd2~=Inf
    [alpha,betta] = circcirc(c1_alpha,c1_beta,rd1,c2_alpha,c2_beta,rd2);
    
    % Estimation of the metric transformation
    A=[1/abs(betta(1)) -alpha(1)/abs(betta(1)) 0;0 1 0;0 0 1];
    
    % Similarity transformation
    l1S=l1/(A*P);l1S=l1S/l1S(3);
    angle=atan(l1S(1)/l1S(2));

    R = [cos(angle) -sin(angle); sin(angle) cos(angle)];
    S = [R zeros(2,1);0 0 1];
    
    % Metrically rectified image
    H=[4 0 0;0 4 0; 0 0 1]*S*A*P;
else
    H=eye(3);
    s=-1;
    disp('Warning: Homography = I');
end

p1_=H*p1';p1_=p1_/p1_(3);
p2_=H*p2';p2_=p2_/p2_(3);
p3_=H*p3';p3_=p3_/p3_(3);
p4_=H*p4';p4_=p4_/p4_(3);
R = makeresampler({'cubic','nearest'},'bound');
Im_metr = imtransform(im,(maketform('projective',(H'))),R,'FillValues',[10;20;30],'XYScale',1);

% show metric result
imshow(Im_metr);
