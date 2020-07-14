%%
ieee='b';
accuracy='real*8';
%% topography
%% X = 400km Y = 300km (dx = dy = 4km)
xx = 1:1:100;
yy = 1:1:75;
[ym,xm] = meshgrid(yy,xx);
x0 = 50; y0 = 37.5;

Ls = 6; % scale of seamount (24 km with dx=4 km)

rsq = (xm-x0).^2+(ym-y0).^2; % r^2
dist0 = exp(-rsq/Ls.^2);
H = 4500;
hm = 0.9*4500;
Ho = hm*dist0-H;

Ho(:,1) = 0;
Ho(:,end) = 0;

fid=fopen('topog_seam1c_seam_clo2.dat','w',ieee); fwrite(fid,Ho,accuracy); fclose(fid);

H0 = Ho*0-4500;
H0(:,1) = 0;
H0(:,end) = 0;
fid=fopen('topog_seam1c_flat_clo2.dat','w',ieee); fwrite(fid,H0,accuracy); fclose(fid);
% % Wind-stress
% tauMax=0.1;
% x=((1:nx)-0.5)/(nx-1); % nx-1 accounts for a solid wall
% y=((1:ny)-0.5)/(ny-1); % ny-1 accounts for a solid wall
% [X,Y]=ndgrid(x,y);
% tau=tauMax*sin(pi*Y);
% fid=fopen('windx.sin_y','w',ieee); fwrite(fid,tau,accuracy); fclose(fid);

%% bnd condition
dy = 4000;
dx = 4000;
x = xx*dx;
y = yy*dy;
z = 0:-50:-500;
z_deep = [z,-1000:-500:-4500];

um1 = 0.1; um2 = 0.1;
hs = 120;  drho = 2;
hc = 200;  ht = 100;
hd = 80;

f = 2*7.29e-5*sind(17);
g = 9.81;
rho0 = 1020;

u1d = um1 + um2*tanh((z+hs)/hd);
u1d_deep = [u1d,0,0,0,0,0,0,0];
for iy = 1:length(y)
    u3d(iy,:) = 0.5*(u1d_deep(1:end-1)+u1d_deep(2:end));
end

rho1d = rho0 - drho*tanh((z+hc)/ht);
dudz = (u1d(1:end-1)-u1d(2:end))/50;
for iy = 1:length(y)
    rho_therm(iy,:) = rho0/g*f*dudz*(y(iy)-y(1));
    rho3d(iy,:) = 0.5*(rho1d(1:end-1)+rho1d(2:end))+rho_therm(iy,:);
end

rho_ref = mean(rho3d,1);
rho_loc = rho3d-rho_ref;
rhoNil = 998.0;
tAlpha = 2e-4;
tref = 20;

t1d = tref+drho/tAlpha/rhoNil*tanh((z+hc)/ht);
for iy = 1:length(y)
    t_therm(iy,:) = -rho0/tAlpha/rhoNil/g*f*dudz*(y(iy)-y(1));
    t3d(iy,:) = 0.5*(t1d(1:end-1)+t1d(2:end))+t_therm(iy,:);
end
t3d_lower(:,1) = t3d(:,end)-0.01;
for iz = 2:8
    t3d_lower(:,iz) = t3d_lower(:,iz-1)-0.1;
end
t3d_deep = cat(2,t3d,t3d_lower);

temp_ref = mean(t3d_deep,1);
tmer1 = t3d_deep(1,:);
tmer2 = t3d_deep(end,:);
for ix = 1:length(x)
    tmer_north(ix,:) = tmer2;
    tmer_south(ix,:) = tmer1;
    vmer_3d(ix,:) = 0*tmer1;
end

% u1 = 0.25;
% dy = 400;
% y = dy*y0; 
% eta = -f/g*(y-y(1))*u1;
% eta = eta';
%%
% fid=fopen('OB_WestH.bin','w',ieee); fwrite(fid,eta,accuracy); fclose(fid);
 fid=fopen('OBzonalU_velo2.bin','w',ieee); fwrite(fid,u3d,accuracy); fclose(fid);
 fid=fopen('OBzonalT_velo2.bin','w',ieee); fwrite(fid,t3d_deep,accuracy); fclose(fid);
%  fid=fopen('OB_northT_velo2.bin','w',ieee); fwrite(fid,tmer_north,accuracy); fclose(fid);
%  fid=fopen('OB_southT_velo2.bin','w',ieee); fwrite(fid,tmer_south,accuracy); fclose(fid);
%  fid=fopen('OB_meriV_velo2.bin','w',ieee); fwrite(fid,vmer_3d,accuracy); fclose(fid);
