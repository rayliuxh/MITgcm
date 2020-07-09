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

spd_p = sqrt(up.^2+vp.^2);

state2_file = '../seamount_0f_flat/state.0000000000.t001.nc';
eta2 = ncread(state2_file,'Eta');
t2 = ncread(state2_file,'Temp');
u2 = ncread(state2_file,'U');
v2 = ncread(state2_file,'V');
size(eta2,3)
up2 = 0.5*(u2(1:end-1,:,:,:)+u2(2:end,:,:,:));
vp2 = 0.5*(v2(:,1:end-1,:,:)+v2(:,2:end,:,:));
spd_p2 = sqrt(up2.^2+vp2.^2);

%% plotting
vscale = 0.3;
vskp = 6;
figure, 
for it = 1:5:size(eta,3)
contourf(x,y,t(:,:,1,it)',28:0.25:31,'linestyle','none')
% contourf(x,y,eta(:,:,it)','linestyle','none')
set(gca,'clim',[28 30])
hold on
quiver(x(1:vskp:end),y(1:vskp:end),...
    up(1:vskp:end,1:vskp:end,1,it)'*vscale,vp(1:vskp:end,1:vskp:end,1,it)'*vscale,0,'k')
title(['Temperature at t=',num2str(it)])
colorbar
pause
end
%% plotting
t(t==0) = nan;
up(up==0) = nan;
vp(vp==0) = nan;
vscale = 10;
vskp = 3;
figure, 
for it = 1:1:size(eta,3)
contourf(x,y,t(:,:,11,it)','linestyle','none')
% pcolor(x(1:vskp:end),y(1:vskp:end),t(1:vskp:end,1:vskp:end,11,it)')
% contourf(x,y,eta(:,:,it)','linestyle','none')
hold on
quiver(x(1:vskp:end),y(1:vskp:end),up(1:vskp:end,1:vskp:end,11,it)'*vscale,...
    vp(1:vskp:end,1:vskp:end,11,it)'*vscale,0,'k')
title(['Temperature at t=',num2str(it)])
colorbar
pause
end

%%
itmax = size(up2,4);
dup = up(:,:,:,itmax)-up2;
dvp = vp(:,:,:,itmax)-vp2;
dspd_p = spd_p(:,:,:,itmax)-spd_p2;
dup(up(:,:,:,1:itmax)==0) = nan;
dvp(vp(:,:,:,1:itmax)==0) = nan;
dspd_p(spd_p(:,:,:,1:itmax)==0) = nan;
spd_p(spd_p==0)=nan;

ik = 12;
vscale = 100;
vskp = 3;
figure, 
% for it = 1:30:size(eta,3)
it = 26;
contourf(x,y,spd_p(:,:,ik,it))
hold on
quiver(x(1:vskp:end),y(1:vskp:end),up(1:vskp:end,1:vskp:end,ik,it)'*vscale,...
    vp(1:vskp:end,1:vskp:end,ik,it)'*vscale,0,'k')
title(['Speed at t=',num2str(it)])
colorbar
% pause
% end
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
