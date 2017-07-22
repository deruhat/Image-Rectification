%% Author: Abdulellah Abualshour
%  King Abdullah University of Science and Technology

clear all
close all

A = imread('ggg.jpg');
I = rgb2gray(A);

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

% display image with lines detected
figure, imshow(I), hold on
for k = 1:size(L,1)
   xy = [L(k,1:2);L(k,3:4)];
   plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
end

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

% visualization
pred_y = polyval(best_model,best_inliers(1,:));
figure;
plot(x,y, 'ro');
hold on; 
plot(best_inliers(1,:),pred_y, '-', 'LineWidth',2);
title('Line fitting found with RANSAC');
legend('Data', 'Model');


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
newImage = imtransform(A, tform);

figure;
imshow(newImage);

A=L(:,1:2)';B=L(:,3:4)';k=size(A,2);
A_=Hv*[A;ones(1,k)];A_=A_./repmat(A_(3,:),[3 1]);
B_=Hv*[B;ones(1,k)];B_=B_./repmat(B_(3,:),[3 1]);
Lv=[A_(1,:)-x;A_(2,:);B_(1,:)-x;B_(2,:)]';

figure, imshow(newImage), hold on
for k = 1:size(L,1)
   xy = [Lv(k,1:2);Lv(k,3:4)];
   plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
end

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

% display image with lines detected
figure, imshow(newImage), hold on
for k = 1:size(Lv,1)
   xy = [Lv(k,1:2);Lv(k,3:4)];
   plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
end

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

% visualization
pred_y_h = polyval(best_model,best_inliers_h(1,:));
figure;
plot(x,y, 'ro');
hold on; 
plot(best_inliers_h(1,:),pred_y_h, '-', 'LineWidth',2);
title('Line fitting found with RANSAC');
legend('Data', 'Model');

figure, imshow(newImage), hold on
max_len = 0;
for k = 1:size(inliers_h,2)
   xy = [inliers_x1_h(k) inliers_y1_h(k);inliers_x2_h(k) inliers_y2_h(k)];
   plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
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

figure;
imshow(finalImage);






