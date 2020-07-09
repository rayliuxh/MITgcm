%%
ieee='b';
accuracy='real*8';
%% topography
x0 = 1:1:40;
y0 = 1:1:40;
[xm0,ym0] = meshgrid(x0,y0);

% Ls = 25; %# 20 km with dx=0.4 km
% dist0 = exp(-((xm0-120).^2+(ym0-120).^2)/Ls.^2);
% Ho = 3700*dist0-4000;
Ho = xm0*0-4000;

fid=fopen('topog_seam0c.dat','w',ieee); fwrite(fid,Ho,accuracy); fclose(fid);

% % Wind-stress
% tauMax=0.1;
% x=((1:nx)-0.5)/(nx-1); % nx-1 accounts for a solid wall
% y=((1:ny)-0.5)/(ny-1); % ny-1 accounts for a solid wall
% [X,Y]=ndgrid(x,y);
% tau=tauMax*sin(pi*Y);
% fid=fopen('windx.sin_y','w',ieee); fwrite(fid,tau,accuracy); fclose(fid);

%% bnd condition
f = 2*7.29e-5*sind(17);
g = 9.81;

u1 = 0.25;

dy = 400;
y = dy*y0;

eta = -f/g*(y-y(1))*u1;
eta = eta';

% fid=fopen('OB_WestH.bin','w',ieee); fwrite(fid,eta,accuracy); fclose(fid);
