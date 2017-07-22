%% Author: Abdulellah Abualshour
%  King Abdullah University of Science and Technology

function [f,H] = rect_func_v(d,length,width,A,B)

d1 = d(1);
d2 = d(2);
w = width;
l = length;

H = [1+((d2-d1)/w)    -d1/l       d1;
     0            1+((d2-d1)/w)   0;
     0            ((d2-d1)/(w*l)) 1];

s=max(1,w/(w+d(2)-d(1)));

S = [s 0 0;0 s 0;0 0 1];

H=S*H;

if isempty(A) && isempty(B)
    
    f=0;
    
else
    
    k=size(A,2);
    
    A_=H*[A;ones(1,k)];A_=A_./repmat(A_(3,:),[3 1]);
    
    B_=H*[B;ones(1,k)];B_=B_./repmat(B_(3,:),[3 1]);
    
    D_=B_(1:2,:)-A_(1:2,:);
    
    f=sum(abs(D_(1,:)));
end

end