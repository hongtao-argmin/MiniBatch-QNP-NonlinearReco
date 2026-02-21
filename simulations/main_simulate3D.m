%% Scattered field from a bead using Lippmann-Schwinger model
clear;
close all;
useGPU(0);
gpuDevice(1);
%% Create two beads
simP = struct();
simP.realData = false;
%cd ~/Tao
simP.ObjectType = 'sheppBorn';%'shepp'; 'Tao'
strM = 'lipp';%'Born';%
simP.lambda0 = 0.406;% wavelength
simP.Nx = 64;%Number of voxels X
simP.Ny = 64;% Number of voxels Y
simP.Nz = 64;% Number of voxels Z
simP.Nroi = [simP.Nx,simP.Ny,simP.Nz];

simP.Lx = simP.Nroi(1)*0.05;%0.05;% Physical length X
simP.Ly = simP.Nroi(2)*0.05;%*0.05;% Physical length Y
simP.Lz = simP.Nroi(3)*0.05;%*0.05;% Physical length Z
simP.n0 = 1.333;% Background refractive index

%radii = [simP.lambda0,simP.lambda0];% Beads radii
%centers = cat(1,[simP.Nx,simP.Ny,simP.Nz/2]/2,...
%    [simP.Nx,simP.Ny,simP.Nz*3/2]/2);% Beads center
%[0.08];% 0.05 refractive index (diff) of the beads
if strcmpi(strM,'Born')
    simP.dn = 0.03;
elseif strcmpi(strM,'lipp')
    %simP.dn = 0.08;
    %simP.dn = 0.1;
    simP.dn = 0.09;
    %simP.dn = 0.11;
    %simP.dn = 0.12;
    %simP.dn = 0.13;
    %simP.dn = 0.14;
    %simP.dn = 0.15;
    %simP.dn = 0.16;
    %simP.dn = 0.17;
    %simP.dn = 0.18;
    %simP.dn = 0.19;
    %simP.dn = 0.2;
end

simP.dx = simP.Lx/simP.Nx;% step size X
simP.dy = simP.Ly/simP.Ny;% step size Y
simP.dz = simP.Lz/simP.Nz;% step size Z
simP.k0 = 2*pi/simP.lambda0;% Freespace wavenumber
simP.k = simP.k0*simP.n0;% Wavenumber
assert(simP.dx==simP.dy & simP.dy==simP.dz,'Step sizes must be equal');

simP.NA = simP.n0;
simP.Ltheta = {deg2rad([-42,42]),deg2rad([-42,42])};
simP.Lphantom = [simP.Lx,simP.Ly,simP.Lz];
if strcmpi(simP.ObjectType,'sheppBorn')
    n = getSheppLogan3D4Born(simP);
elseif strcmpi(simP.ObjectType,'shepp')
    n = getSheppLogan3D(simP);
elseif strcmpi(simP.ObjectType,'Tao')
    n  = ThreeDPhantonImage(simP.Nx,simP.Ny,simP.Nz,simP.n0,simP.dn);
end

f = gpuCpuConverter((simP.k*simP.dx)^2*((n/simP.n0).^2 - 1));% Scattering Potential

simP.Ntheta = 60;

simP.thetas = max(simP.Ltheta{1})*[cos(2*pi*(1:simP.Ntheta)'/simP.Ntheta),...
    sin(2*pi*(1:simP.Ntheta)'/simP.Ntheta)];%Angle of incidence (rad.)

SP.Nmeas = [512,512];%Number of sensors

%Forward model
simP.Niter = 1000;% # iterations to compute the forward model
simP.xtol = 1e-9;% tolerance for convergence criteria
%% Simulate measurements for all angles (ideal plane wave)
Vectx = (1:simP.Nx) - simP.Nx/2;
Vecty = (1:simP.Ny) - simP.Ny/2;
Vectz = (1:simP.Nz) - simP.Nz/2;
[X,Y,Z] = meshgrid(Vectx*simP.dx,Vecty*simP.dy,Vectz*simP.dz);

%XY positions: origin at the center of the sample
[SP.x,SP.y] = meshgrid((simP.dx:simP.dx:simP.dx*SP.Nmeas(1)) - simP.dx*SP.Nmeas(1)/2,...
    (simP.dy:simP.dy:simP.dy*SP.Nmeas(2)) - simP.dy*SP.Nmeas(2)/2);

SP.y = SP.y(:); SP.x = SP.x(:);

%one pixel outside of the sample.
SP.z = simP.dz*(simP.Nroi(3)/2+1)*ones(length(SP.y),1);
simP.Nx = SP.Nmeas(1);
simP.Ny = SP.Nmeas(2);
simP.Lx = simP.Nx*simP.dx;
simP.Ly = simP.Ny*simP.dy;
%% Build the forward model based on Lippmann-Schwinger

% Green function in the volume
G = FreeSpaceKernel3D(simP.k*simP.dx,simP.Nroi,'Vainikko');
% Green function to get the scattered field at the sensor positions
Gtil = Vainikkotilde(simP.k,simP.Nroi,SP,[simP.dx,simP.dy,simP.dz],ones(1,3));

if strcmpi(strM,'lipp')
    H = OpLipp3D(G,Gtil,LinOpIdentity(simP.Nroi),zeros(simP.Nroi),simP.Niter,simP.xtol,[],false);
elseif strcmpi(strM,'born')
    H = LinOpBorn3D(Gtil,LinOpIdentity(simP.Nroi),zeros(simP.Nroi),[],false);
end
% (Optional) Repropagate the scattered field at the (center + shiftZ)
par.distuM = -SP.z(end);%-(simP.Nz)/2*simP.dz;%par.distUin = par.distuM;

Vectx = (1:SP.Nmeas(1)) - SP.Nmeas(1)/2;
Vecty = (1:SP.Nmeas(2)) - SP.Nmeas(2)/2;
[Xmeas,Ymeas] = meshgrid(Vectx*simP.dx,Vecty*simP.dy);

%incident field at the repropagated z-position.
%Propagator
P = LinOpTiltedPropagator(simP.lambda0,simP.n0,par.distuM,simP.dx,SP.Nmeas,[0,0],'AS');
uM.GTc = zeros([H.sizeout,simP.Ntheta]);
uem = zeros([H.sizeout,simP.Ntheta]);
%%
for kk = 1:simP.Ntheta    
    
    simP.kvec = simP.k*[sin(simP.thetas(kk,1));sin(simP.thetas(kk,2))];%Wave vector
    simP.kvec(3) = (-1)^((abs(simP.thetas(kk,1)) > pi/2) + (abs(simP.thetas(kk,2)) > pi/2))*sqrt(simP.k^2 - sum(simP.kvec(1:2).^2));
    uin = gpuCpuConverter(exp(1i*(simP.kvec(1)*X + simP.kvec(2)*Y + simP.kvec(3)*Z)));
    H.uin = uin;
    % Compute the scattered field on the sensors positions
    tic
    y = H*f;
    toc
    % Orthoviews of 3D total field
    %myOrthoviews(abs(H.u));

    yprop = P*y;
    %%
    cuem = exp(1i*(simP.kvec(1)*Xmeas + simP.kvec(2)*Ymeas + 0*simP.kvec(3)*simP.dz));
    if 0
        myOrthoviews(abs(H.u));
        figure(1);
        subplot(221);
        imagesc(abs(y));axis image;colorbar;title('scattered field');
        subplot(222);
        imagesc(abs(yprop));axis image;colorbar;title('centered scattered field');
        subplot(223);
        imagesc(abs(yprop+cuem));axis image;colorbar;title('centered total field');
        subplot(224);
        imagesc(angle((yprop+cuem)./cuem));axis image;colorbar;title('centered relative phase');
        drawnow;
    end
    %%
    uM.GTc(:,:,kk) = gather(yprop + cuem);
    uem(:,:,kk) = gather(cuem);
    fprintf('Approximate expected relative phase %1.2f\n',...
        simP.dx*simP.k0*max(max(sum(n-simP.n0,3))));
end
save(sprintf('%s_%d_HalfAnglesimulation_%s_model_CG_n0_%1.2f_dn_%1.2f_Nviews_%i_Nx_%i_Ny_%i_Nz_%i_dx_%1.2f_dy_%1.2f_dz_%1.2f_GaussianBeam_0_photon_1_illum_circle_%i_%i.mat',...
    strM,simP.Ntheta,simP.ObjectType,simP.n0,simP.dn(1),simP.Ntheta,simP.Nx,simP.Ny,simP.Nz,simP.dx,simP.dy,simP.dz,rad2deg(simP.Ltheta{1})),...
    'SP','n','simP','uM','uem','-v7.3');
