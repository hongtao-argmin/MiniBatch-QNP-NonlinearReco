function [xIndexSet,yIndexSet]= draw_rectangle(center_location,L,H,theta,N)
% Inputs: 
%   center_location: center location of the ractangle.
%   L,H: the length.
%   theta: rotate angle. e.g., pi/4.
%   N: the number of sampling points.
% Outputs: xIndexSet,yIndexSet: the obtained four indeces.
% 
% the domain is set to be [-1 1];
center_location = center_location - center_location/2;
center1 = center_location(1);
center2 = center_location(2);
R = ([cos(theta), -sin(theta); sin(theta), cos(theta)]); 
X = ([-L/2, L/2, L/2, -L/2]);
Y = ([-H/2, -H/2, H/2, H/2]);

for i=1:4
    T(:,i)=R*[X(i); Y(i)];
end

x_lower_left = center1+T(1,1);
x_lower_right = center1+T(1,2);
x_upper_right = center1+T(1,3);
x_upper_left = center1+T(1,4);

y_lower_left = center2+T(2,1);
y_lower_right = center2+T(2,2);
y_upper_right = center2+T(2,3);
y_upper_left = center2+T(2,4);

x_coor = [x_lower_left x_lower_right x_upper_right x_upper_left];
y_coor = [y_lower_left y_lower_right y_upper_right y_upper_left];

xIndexSet = round(x_coor.*(N/2))+N/2;
yIndexSet = round(y_coor.*(N/2))+N/2;

index = xIndexSet<0;
xIndexSet(index)=0;
index = xIndexSet>N;
xIndexSet(index) = N;

index = yIndexSet<0;
yIndexSet(index)=0;
index = yIndexSet>N;
yIndexSet(index) = N;

return
