% Generate synthetic 3D image
%
% Author: Tao  Hong 8th Dec.
function n = ThreeDPhantonImage(M,N,K,n_0,dn)
% M = 64;
% N = 64;
% K = 64;
% n_0 = 1.33;
% dn = 0.1;

n = n_0*ones(M,N,K);
partial = 0.03;
M_begin = round(M*partial);
M_end = M-round(M*partial);
N_begin = round(N*partial);
N_end = N-round(N*partial);
K_begin = round(K*partial);
K_end = K-round(K*partial);


n(M_begin:M_end,N_begin:N_end,K_begin:K_end) = n(M_begin:M_end,N_begin:N_end,K_begin:K_end)+dn;
% add circle along the z-direction with increase and decrease radius.
% add big out circle
center_pos_1 = 0.5;
center_pos_2 = 0.5;

center_1 = round(center_pos_1*M);
center_2 = round(center_pos_2*N);
[X,Y] = meshgrid(1:M,1:N);
radius_1 = 0.42;
radius_2 = 0.08;
begin_index = 5;
dn_radius = (radius_1-radius_2)/(round(K/2)-begin_index);
for iter = begin_index:K-begin_index
    if iter<=round(K/2)
        radius_N = round((radius_2+dn_radius*(iter-begin_index))*N);
    else
        radius_N = round((radius_1-dn_radius*(iter-round(K/2)))*N);
    end
    temp =  squeeze(n(:,:,iter));
    temp((X-center_1).^2+(Y-center_2).^2<=radius_N^2) =  temp((X-center_1).^2+(Y-center_2).^2<=radius_N^2)-0.2*dn;
    n(:,:,iter) = temp;
end

% put circle in Y direction
radius = 0.2;
radius_N = round(radius*N);
for iter = begin_index:N-begin_index
    temp =  squeeze(n(:,iter,:));
    temp((X-center_1).^2+(Y-center_2).^2<=radius_N^2) =  temp((X-center_1).^2+(Y-center_2).^2<=radius_N^2)-0.3*dn;
    n(:,iter,:) = temp;
end

% put four shape in X direction. equal size
for iter = begin_index:M-begin_index
    temp =  squeeze(n(iter,:,:));
    [xIndexSet,yIndexSet]= draw_rectangle([0 0],0.35,0.35,0,N);
    mask = poly2mask(xIndexSet, yIndexSet, M, M);
    temp(mask) =  temp(mask)-0.35*dn;
    [xIndexSet,yIndexSet]= draw_rectangle([-1.5 -1.5],0.25,0.35,0,N);
    mask = poly2mask(xIndexSet, yIndexSet, M, M);
    temp(mask) =  temp(mask)-0.3*dn;

    [xIndexSet,yIndexSet]= draw_rectangle([1.5 1.5],0.25,0.25,0,N);
    mask = poly2mask(xIndexSet, yIndexSet, M, M);
    temp(mask) =  temp(mask)-0.3*dn;
    n(iter,:,:) = temp;
end


% strategy III
% [xIndexSet,yIndexSet]= draw_rectangle([-0.93 0.93],1.8,0.2,pi/4,N);
% mask = poly2mask(xIndexSet, yIndexSet, N, N);
% mask = mask&(X-center_1).^2+(Y-center_2).^2<=radius_N_1^2;
% n(mask) = n(mask)-0.2*dn;
%
% [xIndexSet,yIndexSet]= draw_rectangle([0.93 -0.93],1.8,0.2,pi/4,N);
% mask = poly2mask(xIndexSet, yIndexSet, N, N);
% mask = mask&(X-center_1).^2+(Y-center_2).^2<=radius_N_1^2;
% n(mask) = n(mask)-0.2*dn;
%
% [xIndexSet,yIndexSet]= draw_rectangle([0 0],2,0.12,-pi/4,N);
% mask = poly2mask(xIndexSet, yIndexSet, N, N);
% mask = mask&(X-center_1).^2+(Y-center_2).^2<=radius_N_1^2;
% n(mask) = n(mask)-0.3*dn;
%
% % plot two circles and two inside circles
% radius_2 = 0.13;
% radius_2_inside = 0.088;
% center_pos_3 = 0.3;
% center_pos_4 = 0.3;
% radius_N_2 = round(radius_2*N);
% radius_N_2_inside = round(radius_2_inside*N);
% center_3 = round(center_pos_3*N);
% center_4 = round(center_pos_4*N);
% n((X-center_3).^2+(Y-center_4).^2<=radius_N_2^2) = n((X-center_3).^2+(Y-center_4).^2<=radius_N_2^2)-0.2*dn;
% n((X-center_3).^2+(Y-center_4).^2<=radius_N_2_inside^2) = n((X-center_3).^2+(Y-center_4).^2<=radius_N_2_inside^2)-0.2*dn;
%
%
% radius_2 = 0.13;
% radius_2_inside = 0.05;
% center_pos_3 = 0.7;
% center_pos_4 = 0.7;
% radius_N_2 = round(radius_2*N);
% radius_N_2_inside = round(radius_2_inside*N);
% center_3 = round(center_pos_3*N);
% center_4 = round(center_pos_4*N);
% n((X-center_3).^2+(Y-center_4).^2<=radius_N_2^2) = n((X-center_3).^2+(Y-center_4).^2<=radius_N_2^2)-0.2*dn;
% n((X-center_3).^2+(Y-center_4).^2<=radius_N_2_inside^2) = n((X-center_3).^2+(Y-center_4).^2<=radius_N_2_inside^2)-0.2*dn;
%


%{
% strategy I 
% plot rotated rectangle
% negative positive: left down corner. x(left): negative. y(down): positive
[xIndexSet,yIndexSet]= draw_rectangle([-0.93 0.93],1.8,0.17,pi/4,N);
mask = poly2mask(xIndexSet, yIndexSet, N, N);
mask = mask&(X-center_1).^2+(Y-center_2).^2<=radius_N_1^2;
n(mask) = n(mask)-0.2*dn;

[xIndexSet,yIndexSet]= draw_rectangle([0.93 -0.93],1.8,0.17,pi/4,N);
mask = poly2mask(xIndexSet, yIndexSet, N, N);
mask = mask&(X-center_1).^2+(Y-center_2).^2<=radius_N_1^2;
n(mask) = n(mask)-0.2*dn;

[xIndexSet,yIndexSet]= draw_rectangle([0 0],2,0.13,pi/4,N);
mask = poly2mask(xIndexSet, yIndexSet, N, N);
mask = mask&(X-center_1).^2+(Y-center_2).^2<=radius_N_1^2;
n(mask) = n(mask)-0.3*dn;

[xIndexSet,yIndexSet]= draw_rectangle([0 0],2,0.12,-pi/4,N);
mask = poly2mask(xIndexSet, yIndexSet, N, N);
mask = mask&(X-center_1).^2+(Y-center_2).^2<=radius_N_1^2;
n(mask) = n(mask)-0.3*dn;


[xIndexSet,yIndexSet]= draw_rectangle([0.7 0.7],1.08,0.12,-pi/4,N);
mask = poly2mask(xIndexSet, yIndexSet, N, N);
mask = mask&(X-center_1).^2+(Y-center_2).^2<=radius_N_1^2;
n(mask) = n(mask)-0.3*dn;


[xIndexSet,yIndexSet]= draw_rectangle([-0.7 -0.7],1.08,0.13,-pi/4,N);
mask = poly2mask(xIndexSet, yIndexSet, N, N);
mask = mask&(X-center_1).^2+(Y-center_2).^2<=radius_N_1^2;
n(mask) = n(mask)-0.3*dn;

% [xIndexSet,yIndexSet]= draw_rectangle([0.5 0.4],2,0.12,-pi/4,N);
% mask = poly2mask(xIndexSet, yIndexSet, N, N);
% mask = mask&(X-center_1).^2+(Y-center_2).^2<=radius_N_1^2;
% n(mask) = n(mask)-0.3*dn;
% [xIndexSet,yIndexSet]= draw_rectangle([-0.3 -0.5],2.3,0.16,-pi/4,N);
% mask = poly2mask(xIndexSet, yIndexSet, N, N);
% mask = mask&(X-center_1).^2+(Y-center_2).^2<=radius_N_1^2;
% n(mask) = n(mask)-0.3*dn;
%}

% strategy II
% add two circles
% radius_2 = 0.05;
% center_pos_3 = 0.3;
% center_pos_4 = 0.3;
%
% radius_N_2 = round(radius_2*N);
% center_3 = round(center_pos_3*N);
% center_4 = round(center_pos_4*N);
% n((X-center_3).^2+(Y-center_4).^2<=radius_N_2^2) = n((X-center_3).^2+(Y-center_4).^2<=radius_N_2^2)-0.3*dn;


% radius_2 = 0.06;
% center_pos_3 = 0.7;
% center_pos_4 = 0.7;
%
% radius_N_2 = round(radius_2*N);
% center_3 = round(center_pos_3*N);
% center_4 = round(center_pos_4*N);
% n((X-center_3).^2+(Y-center_4).^2<=radius_N_2^2) = n((X-center_3).^2+(Y-center_4).^2<=radius_N_2^2)-0.25*dn;



%{
% add first small circle
% radius_2 = 0.045;
% center_pos_3 = 0.2;
% center_pos_4 = 0.19;
% 
% radius_N_2 = round(radius_2*N);
% center_3 = round(center_pos_3*N);
% center_4 = round(center_pos_4*N);
% n((X-center_3).^2+(Y-center_4).^2<=radius_N_2^2) = n((X-center_3).^2+(Y-center_4).^2<=radius_N_2^2)-0.6*dn;

% % add first inside circle
% radius_3 = 0.05;
% radius_N_3 = round(radius_3*N);
% n((X-center_3).^2+(Y-center_4).^2<=radius_N_3^2) = n((X-center_3).^2+(Y-center_4).^2<=radius_N_3^2)+0.1*dn;


% add one out rectangle
% x_corrdinate = round([760 1000 62 965]/1024*N); % (y,y,x,x)
% % index = (Y>=x_corrdinate(1)&Y<=x_corrdinate(2)&X>=x_corrdinate(3)&X<=x_corrdinate(4));
% % n(index) = n(index)+dn;
% index = (Y>=x_corrdinate(1)&Y<=x_corrdinate(2)&X>=x_corrdinate(3)&X<=x_corrdinate(4)&((X-center_1).^2+(Y-center_2).^2>radius_N_1^2));
% n(index) = n(index)+0.8*dn;

% add inside rectangle
x_corrdinate = round([750 850 380 680]/1024*N); % (y,y,x,x)
index = (Y>=x_corrdinate(1)&Y<=x_corrdinate(2)&X>=x_corrdinate(3)&X<=x_corrdinate(4));
n(index) = n(index)-0.25*dn;

% add several small circles but with major difference
radius_1 = 0.05;
centerSmall_pos_1 = 0.89;
centerSmall_pos_2 = 0.35;
radiusSmall_N_1 = round(radius_1*N);
centerSmall_1 = round(centerSmall_pos_1*N);
centerSmall_2 = round(centerSmall_pos_2*N);
n((X-centerSmall_1).^2+(Y-centerSmall_2).^2<=radiusSmall_N_1^2) = n((X-centerSmall_1).^2+(Y-centerSmall_2).^2<=radiusSmall_N_1^2)-0.4*dn;

radius_1 = 0.06;
centerSmall_pos_1 = 0.88;
centerSmall_pos_2 = 0.63;
radiusSmall_N_1 = round(radius_1*N);
centerSmall_1 = round(centerSmall_pos_1*N);
centerSmall_2 = round(centerSmall_pos_2*N);
n((X-centerSmall_1).^2+(Y-centerSmall_2).^2<=radiusSmall_N_1^2) = n((X-centerSmall_1).^2+(Y-centerSmall_2).^2<=radiusSmall_N_1^2)-0.5*dn;

% radius_1 = 0.09;
% centerSmall_pos_1 = 0.75;
% centerSmall_pos_2 = 0.5;
% radiusSmall_N_1 = round(radius_1*N);
% centerSmall_1 = round(centerSmall_pos_1*N);
% centerSmall_2 = round(centerSmall_pos_2*N);
% n((X-centerSmall_1).^2+(Y-centerSmall_2).^2<=radiusSmall_N_1^2) = n((X-centerSmall_1).^2+(Y-centerSmall_2).^2<=radiusSmall_N_1^2)-0.2*dn;

radius_1 = 0.05;
centerSmall_pos_1 = 0.53;
centerSmall_pos_2 = 0.47;
radiusSmall_N_1 = round(radius_1*N);
centerSmall_1 = round(centerSmall_pos_1*N);
centerSmall_2 = round(centerSmall_pos_2*N);
n((X-centerSmall_1).^2+(Y-centerSmall_2).^2<=radiusSmall_N_1^2) = n((X-centerSmall_1).^2+(Y-centerSmall_2).^2<=radiusSmall_N_1^2)-0.5*dn;

radius_1 = 0.06;
centerSmall_pos_1 = 0.66;
centerSmall_pos_2 = 0.145;
radiusSmall_N_1 = round(radius_1*N);
centerSmall_1 = round(centerSmall_pos_1*N);
centerSmall_2 = round(centerSmall_pos_2*N);
n((X-centerSmall_1).^2+(Y-centerSmall_2).^2<=radiusSmall_N_1^2) = n((X-centerSmall_1).^2+(Y-centerSmall_2).^2<=radiusSmall_N_1^2)-0.5*dn;
% add one inside triangle
xCoords = round(([536 268 650]/1024)*N);
yCoords = round(([31 293 255]/1024)*N);
mask = poly2mask(xCoords, yCoords, N, N);
n(mask) = n(mask)-0.5*dn;%2*dn/3;

radius_1 = 0.08;
centerSmall_pos_1 = 0.12;
centerSmall_pos_2 = 0.55;
radiusSmall_N_1 = round(radius_1*N);
centerSmall_1 = round(centerSmall_pos_1*N);
centerSmall_2 = round(centerSmall_pos_2*N);
n((X-centerSmall_1).^2+(Y-centerSmall_2).^2<=radiusSmall_N_1^2) = n((X-centerSmall_1).^2+(Y-centerSmall_2).^2<=radiusSmall_N_1^2)-0.3*dn;
%}

% % add one inside triangle
% xCoords = round(([740 690 866]/1024)*N);
% yCoords = round(([305 427 429]/1024)*N);
% mask = poly2mask(xCoords, yCoords, N, N);
% n(mask) = n(mask)-0.3*dn;%2*dn/3;0.75*


% % add two larger disks.
%
% radius = 0.25;%0.18;
% radius_1 = 0.13;
% radius_2 = 0.13;
%
% center_pos_1 = 0.3;
% center_pos_2 = 0.5;
%
% center_pos_3 = 0.68;
% center_pos_4 = 0.5;
% radius_N_1 = round(radius*N);
% radius_N_2 = round(radius*N);
% radius_N_3 = round(radius_1*N);
% radius_N_4 = round(radius_2*N);
% center_1 = round(center_pos_1*N);
% center_2 = round(center_pos_2*N);
% center_3 = round(center_pos_3*N);
% center_4 = round(center_pos_4*N);
%
% [X,Y] = meshgrid(1:N,1:N);
% n((X-center_1).^2+(Y-center_2).^2<=radius_N_1^2) = n((X-center_1).^2+(Y-center_2).^2<=radius_N_1^2)+0.8*dn;
%
% n((X-center_3).^2+(Y-center_4).^2<=radius_N_2^2) = n((X-center_3).^2+(Y-center_4).^2<=radius_N_2^2)+0.8*dn;
%
% n((X-center_3).^2+(Y-center_4).^2<=radius_N_2^2&(X-center_1).^2+(Y-center_2).^2<=radius_N_1^2) = n((X-center_3).^2+(Y-center_4).^2<=radius_N_1^2&(X-center_1).^2+(Y-center_2).^2<=radius_N_1^2)-0.8*dn;
%
% n((X-center_1).^2+(Y-center_2).^2<=radius_N_3^2) = n((X-center_1).^2+(Y-center_2).^2<=radius_N_3^2)+0.2*dn;
%
% n((X-center_3).^2+(Y-center_4).^2<=radius_N_4^2) = n((X-center_3).^2+(Y-center_4).^2<=radius_N_4^2)+0.2*dn;
%
%
%
% % add four rectangles
% x_corrdinate = round([40 390 80 480]/1024*N); % (y,y,x,x)
% index = (Y>=x_corrdinate(1)&Y<=x_corrdinate(2)&X>=x_corrdinate(3)&X<=x_corrdinate(4))&((X-center_1).^2+(Y-center_2).^2>radius_N_1^2);
% n(index) = n(index)+0.9*dn;
%
% % x_corrdinate = round([45 390 700 1000]/1024*N); % (y,y,x,x)
% % index = (Y>=x_corrdinate(1)&Y<=x_corrdinate(2)&X>=x_corrdinate(3)&X<=x_corrdinate(4))&((X-center_3).^2+(Y-center_4).^2>radius_N_2^2);
% % n(index) = n(index)+0.85*dn;
% radius_1 = 0.05;%0.18;
% radius_2 = 0.07;
% center_pos_5 = 0.85;
% center_pos_6 = 0.12;
%
% center_pos_7 = 0.13;
% center_pos_8 = 0.88;
%
% radius_N_3 = round(radius_1*N);
% radius_N_4 = round(radius_2*N);
%
% %center_5 = round(center_pos_5*N);
% center_6 = round(center_pos_6*N);
% center_9 = round((center_pos_5-0.22)*N);
%
% center_7 = round(center_pos_7*N);
% center_8 = round(center_pos_8*N);
% center_10 = round((center_pos_7+0.23)*N);
% %n((X-center_5).^2+(Y-center_6).^2<=radius_N_3^2) = n((X-center_5).^2+(Y-center_6).^2<=radius_N_3^2)+dn;
% %n((X-center_9).^2+(Y-center_6).^2<=radius_N_3^2) = n((X-center_9).^2+(Y-center_6).^2<=radius_N_3^2)+dn;
%
% xCoords_1 = round(([820 590 960]/1024)*N);
% yCoords_1 = round(([40 185 185]/1024)*N);
% mask = poly2mask(xCoords_1, yCoords_1, N, N);
% n(mask) = n(mask)+dn;
%
%
% % x_corrdinate = round([630 970 80 280]/1024*N); % (y,y,x,x)
% % index = (Y>=x_corrdinate(1)&Y<=x_corrdinate(2)&X>=x_corrdinate(3)&X<=x_corrdinate(4))&((X-center_1).^2+(Y-center_2).^2>radius_N_1^2);
% % n(index) = n(index)+0.85*dn;
% n((X-center_7).^2+(Y-center_8).^2<=radius_N_4^2) = n((X-center_7).^2+(Y-center_8).^2<=radius_N_4^2)+dn;
% n((X-center_10).^2+(Y-center_8).^2<=radius_N_4^2) = n((X-center_10).^2+(Y-center_8).^2<=radius_N_4^2)+dn;
%
% x_corrdinate = round([630 970 560 900]/1024*N); % (y,y,x,x)
% index = (Y>=x_corrdinate(1)&Y<=x_corrdinate(2)&X>=x_corrdinate(3)&X<=x_corrdinate(4))&((X-center_3).^2+(Y-center_4).^2>radius_N_2^2);
% n(index) = n(index)+0.85*dn;%0.9*dn;%
%
% % xCoords_1 = round(([330 260 420]/1024)*N);
% % yCoords_1 = round(([880 960 950]/1024)*N);
% % mask = poly2mask(xCoords_1, yCoords_1, N, N);
% % n(mask) = n(mask)+dn;
%
% %{
% % plot two symmetric circles
% radius = 0.18;%0.18;
% center_pos = 0.2;
% radius_N = round(radius*N);
% center_1 = round(center_pos*N);
% center_2 = round((center_pos+0.6)*N);
% [X,Y] = meshgrid(1:N,1:N);
% n((X-center_1).^2+(Y-center_1).^2<=radius_N^2) = n((X-center_1).^2+(Y-center_1).^2<=radius_N^2)+dn;
% n((X-center_2).^2+(Y-center_1).^2<=radius_N^2) = n((X-center_2).^2+(Y-center_1).^2<=radius_N^2)+dn;
% % add two more smaller symmetric circles inside these two cycles
% % radius = 0.04;
% % center_pos = 0.25;
% % radius_N = round(radius*N);
% % center_1 = round(center_pos*N);
% % center_2 = round((center_pos+0.5)*N);
% % [X,Y] = meshgrid(1:N,1:N);
% % n((X-center_1).^2+(Y-center_1).^2<=radius_N^2) = n((X-center_1).^2+(Y-center_1).^2<=radius_N^2)+0.35*dn;
% % n((X-center_2).^2+(Y-center_1).^2<=radius_N^2) = n((X-center_2).^2+(Y-center_1).^2<=radius_N^2)+0.35*dn;
%
%
% % add two more smaller symmetric triangleq inside two circles
%
% % xCoords = round(([250 205 300]/1024)*N);
% % yCoords = round(([200 280 280]/1024)*N);
% % mask = poly2mask(xCoords, yCoords, N, N);
% % n(mask) = n(mask)+0.25*dn;
% %
% %
% % xCoords = round(([770 725 822]/1024)*N);
% % yCoords = round(([200 278 278]/1024)*N);
% % mask = poly2mask(xCoords, yCoords, N, N);
% % n(mask) = n(mask)+0.25*dn;
%
%
% % % add triangleq
% % xCoords = round(([255 88 445]/1024)*N);
% % yCoords = round(([155 350 350]/1024)*N);
% % mask = poly2mask(xCoords, yCoords, N, N);
% % n(mask) = n(mask)+dn;%2*dn/3;0.75*
% % % add one rectangle inside the triangleq
% % % x_corrdinate = round([270 305 220 300]/1024*N);
% % % index = Y>=x_corrdinate(1)&Y<=x_corrdinate(2)&X>=x_corrdinate(3)&X<=x_corrdinate(4);
% % % n(index) = n(index)+0.25*dn;
%
% % add small triangleq
% % xCoords = round(([490 310 690]/1024)*N);
% % yCoords = round(([500 600 600]/1024)*N);
% % mask = poly2mask(xCoords, yCoords, N, N);
% % n(mask) = n(mask)+dn;
%
% radius = 0.035;
% center_pos = 0.5;
% radius_N = round(radius*N);
% center_middle_1 = round(center_pos*N);
% center_middle_2 =  round((center_pos-0.2)*N);
% [X,Y] = meshgrid(1:N,1:N);
% n((X-center_middle_1).^2+(Y-center_middle_2).^2<=radius_N^2) = n((X-center_middle_1).^2+(Y-center_middle_2).^2<=radius_N^2)+dn;
%
% radius = 0.085;
% center_pos = 0.5;
% radius_N = round(radius*N);
% center_middle_1 = round(center_pos*N);
% center_middle_2 =  round((center_pos+0.05)*N);
% [X,Y] = meshgrid(1:N,1:N);
% n((X-center_middle_1).^2+(Y-center_middle_2).^2<=radius_N^2) = n((X-center_middle_1).^2+(Y-center_middle_2).^2<=radius_N^2)+dn;
%
% radius = 0.085;
% center_pos = 0.2;
% radius_N = round(radius*N);
% center_middle_1 = round(center_pos*N);
% center_middle_2 =  round((center_pos+0.33)*N);
% [X,Y] = meshgrid(1:N,1:N);
% n((X-center_middle_1).^2+(Y-center_middle_2).^2<=radius_N^2) = n((X-center_middle_1).^2+(Y-center_middle_2).^2<=radius_N^2)+dn;
%
% radius = 0.085;
% center_pos = 0.8;
% radius_N = round(radius*N);
% center_middle_1 = round(center_pos*N);
% center_middle_2 =  round((center_pos-0.25)*N);
% [X,Y] = meshgrid(1:N,1:N);
% n((X-center_middle_1).^2+(Y-center_middle_2).^2<=radius_N^2) = n((X-center_middle_1).^2+(Y-center_middle_2).^2<=radius_N^2)+dn;
%
% % add rectangle
% x_corrdinate = round([760 1000 62 965]/1024*N);
% index = Y>=x_corrdinate(1)&Y<=x_corrdinate(2)&X>=x_corrdinate(3)&X<=x_corrdinate(4);
% n(index) = n(index)+0.75*dn;
% % add two more circles inside the rectangle
% radius = 0.035;
% center_pos_1 = 0.26;
% center_pos_2 = 0.87;
% radius_N = round(radius*N);
% center_1 = round(center_pos_1*N);
% center_2 = round((center_pos_1+0.5)*N);
% center_3 = round(center_pos_2*N);
%
% [X,Y] = meshgrid(1:N,1:N);
% n((X-center_1).^2+(Y-center_3).^2<=radius_N^2) = n((X-center_1).^2+(Y-center_3).^2<=radius_N^2)+0.25*dn;
% n((X-center_2).^2+(Y-center_3).^2<=radius_N^2) = n((X-center_2).^2+(Y-center_3).^2<=radius_N^2)+0.25*dn;
%
% % %add two more separate rectangle
% % x_corrdinate = round([420 500 117 187]/1024*N);
% % index = Y>=x_corrdinate(1)&Y<=x_corrdinate(2)&X>=x_corrdinate(3)&X<=x_corrdinate(4);
% % n(index) = n(index)+0.2*dn;
% %
% % x_corrdinate = round([420 500 327 397]/1024*N);
% % index = Y>=x_corrdinate(1)&Y<=x_corrdinate(2)&X>=x_corrdinate(3)&X<=x_corrdinate(4);
% % n(index) = n(index)+0.2*dn;
% %}