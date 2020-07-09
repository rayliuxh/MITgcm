% This is a matlab script that generates the input data

% Dimensions of grid
nx=40;
ny=40;
nz=9;
% Nominal depth of model (meters)
% Horizontal resolution (m)
dx=30e3;
% Rotation

ieee='b';
prec='real*8';


% initial fields for salinity
ti = [20.,10.,8.,6.,5.5,5.3,5.1,5.0,4.9];
si = [10.,10.,10.,10.,10.,10.,10.,10.,10.];
%fid=fopen('S.init','w',ieee); fwrite(fid,si*ones(nx,ny,nz),prec); fclose(fid);

% open boundary conditions;
u0 = [.25, 0.,0.,0.,0.,0.,0.,0.,0.];

% create two time slabs for testing
for k=1:nz
    uMerid(:,k) = u0(k)*ones(nx,1);
    uZonal(:,k) = u0(k)*ones(ny,1);
    tZonal(:,k) = ti(k)*ones(ny,1);
    sZonal(:,k) = si(k)*ones(ny,1);
end

% to test prescribing Eta in NonLin-FreeSurf formulation:
et1=0.1;
etWest = cat(2,+et1*ones(ny,1),-et1*ones(ny,1));
etEast = cat(2,-et1*ones(ny,1),+et1*ones(ny,1));

%  fid=fopen('OBmeridU.bin','w',ieee); fwrite(fid,uMerid,prec); fclose(fid);
 fid=fopen('OBzonalU.bin','w',ieee); fwrite(fid,uZonal,prec); fclose(fid);
%  fid=fopen('OB_WestU.bin','w',ieee); fwrite(fid,uZonal ,prec); fclose(fid);
%  fid=fopen('OB_EastU.bin','w',ieee); fwrite(fid,uZonal ,prec); fclose(fid);
%  fid=fopen('OBzonalT.bin','w',ieee); fwrite(fid,tZonal,prec); fclose(fid);
%  fid=fopen('OBzonalS.bin','w',ieee); fwrite(fid,sZonal,prec); fclose(fid);
%  fid=fopen('OBzonalW.bin','w',ieee); fwrite(fid,wZonal,prec); fclose(fid);
%  fid=fopen('OB_WestH.bin','w',ieee); fwrite(fid,etWest,prec); fclose(fid);
%  fid=fopen('OB_EastH.bin','w',ieee); fwrite(fid,etEast,prec); fclose(fid);

% %- rbcs mask & restauring tracer field:
% msk=ones(nx,ny,nz);
% xMx=max(x);
% shapeX=(x-xMx)/dx;
% shapeX=exp(shapeX*2/3);
% 
% [I]=find(shapeX < 5.e-3); fprintf('zero out rbc-mask up to i= %i\n',max(I));
% shapeX(I)=0.;
% var=shapeX'*ones(1,ny*nz);
% fid=fopen('rbcs_mask.bin','w',ieee); fwrite(fid,var,prec); fclose(fid);
% 
% tr1=(si+s0)/2;
% var=tr1*ones(nx,ny,nz);
% fid=fopen('rbcs_Tr1_fld.bin','w',ieee); fwrite(fid,var,prec); fclose(fid);

