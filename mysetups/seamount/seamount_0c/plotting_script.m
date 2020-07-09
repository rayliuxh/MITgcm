%%
grd_file = 'grid.t001.nc';
fCori = ncread(grd_file,'fCori');
dxc = ncread(grd_file,'dxC');
%% reading data
state_file = 'state.0000000000.t001.nc';

x = ncread(state_file,'X');
y = ncread(state_file,'Y');
z = ncread(state_file,'Z');
xp = ncread(state_file,'Xp1');
yp = ncread(state_file,'Yp1');

u = ncread(state_file,'U');
v = ncread(state_file,'V');
t = ncread(state_file,'Temp');
eta = ncread(state_file,'Eta');
size(u,4)

up = 0.5*(u(1:end-1,:,:,:)+u(2:end,:,:,:));
vp = 0.5*(v(:,1:end-1,:,:)+v(:,2:end,:,:));

state2_file = 'BACKUP_OB_succ_005_cd07_2ob_kv1e-4/state.0000000000.t001.nc';
eta2 = ncread(state2_file,'Eta');
t2 = ncread(state2_file,'Temp');
size(eta2,3)

%% plotting
vscale = 0.1;
figure, 
for it = 1:30:size(eta,3)
contourf(x,y,t(:,:,1,it)',28:0.25:31,'linestyle','none')
% contourf(x,y,eta(:,:,it)','linestyle','none')
set(gca,'clim',[28 30])
hold on
quiver(x,y,up(:,:,1,it)'*vscale,vp(:,:,1,it)'*vscale,0,'k')
title(['Temperature at t=',num2str(it)])
colorbar
pause
end
%%
figure, 
for it = 1:3:315
   quiver(x,y,up(:,:,1,it)',vp(:,:,1,it)')
   title(['Velocity at t=',num2str(it)])
pause
end
%% for comparison
icom = 600;
figure, 
subplot(121)
contourf(x,y(2:end-1),t(:,2:end-1,1,icom)')
subplot(122)
contourf(x,y(2:end-1),t2(:,2:end-1,1,icom)')
