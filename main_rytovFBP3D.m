%% Main Rytov FBP 3D

par.rytFBP.Lxext = simP.Lx; par.rytFBP.Lyext = simP.Ly;

%if hasandis(simP,'realData')
    %uM.rytFBP = uM.curr_y(:,:,2:);
%else
if ndims(uM.GTc)>3
    uM.rytFBP = squeeze(uM.curr_y(:,:,2,:));
else
    uM.rytFBP = uM.curr_y(:,:,par.curr_ind_thetas);
end


par.rytFBP.roix = round(1 + simP.Nx/2 - simP.Lxroi/simP.dx/2:simP.Nx/2 + simP.Lxroi/simP.dx/2);
par.rytFBP.roiy = round(1 + simP.Ny/2 - simP.Lyroi/simP.dy/2:simP.Ny/2 + simP.Lyroi/simP.dy/2);

if exist('curr_uem','var')
    Muin = curr_uem;
elseif exist('uem','var')
    if ndims(uem)>3
        Muin = squeeze(uem(:,:,2,:));
    else
        Muin = uem;
    end
else
    if ndims(uin) > 3
        Muin = squeeze(uin(:,:,2,:));
    else
        Muin = uin;
    end
end

Muin = Muin(:,:,par.curr_ind_thetas);

if ~simP.realData && ~strcmp(simP.ObjectType,'RBC') && ~contains(simP.ObjectType,'shepp') && ~contains(simP.ObjectType,'printed')
    par.rytFBP.distUin = -curr_SP.z(end);%-simP.Lz/2;
    par.rytFBP.distuM = par.rytFBP.distUin;
    [Muin,uM.rytFBP] = setFields3D(Muin,uM.rytFBP,par.rytFBP,simP,'rytFBP',[],[],false);
    isProp = true;
else
    isProp = false;
end


uM.rytFBP = uM.rytFBP(par.rytFBP.roiy,par.rytFBP.roix,:);



Muin = Muin(par.rytFBP.roiy,par.rytFBP.roix,:);
%%

par.rytFBP.curr_thetas = par.curr_thetas;%true(simP.Ntheta,1);%par.curr_thetas;
par.rytFBP.thetas = simP.thetas(par.rytFBP.curr_thetas,:);

if ~exist('log_imag_ur','var') && ~exist('new_uwr','var') || isfield(par,'rytRef')
    [XX,YY] = meshgrid((1:simP.Nx) - simP.Nx/2, (1:simP.Ny) - simP.Ny/2);
    mask2D = LinOpDiag([simP.Nx,simP.Ny],1);%1*(sqrt(XX.^2 + YY.^2) <= simP.Nx/3));
    if ~isProp
        if isfield(par,'rytRef')
            currR = -par.rytRef;
        else
            currR = 0;
        end
        par.rytFBP.distUin = currR;%-simP.Lz/2;
        par.rytFBP.distuM = par.rytFBP.distUin;
        par.rytFBP.Ntheta = par.Ntheta;
        par.rytFBP.curr_ind_thetas = 1:par.rytFBP.Ntheta;
        [Muin,uM.rytFBP] = setFields3D(Muin,uM.rytFBP,par.rytFBP,simP,'rytFBP',[],[],true);
    end
    uM.rytFBP = squeeze(uM.rytFBP);
    Muin = squeeze(Muin);
    ur = uM.rytFBP./Muin;
    if isfield(par,'rytRref') && exist(sprintf('curr_log_imag_ur_%1.2e.mat',par.rytRef),'file')
        load(sprintf('curr_log_imag_ur_%1.2e.mat',par.rytRef));
    else
        curr_log_imag_ur = zeros(size(ur));
        fprintf('Unwrapping...\n');
        potential.quantized = 'no';
        potential.threshold = pi;
        
        %figure(99);
        for kk = 1:size(uM.rytFBP,3)
            if exist('puma_ho','file')
                curr_log_imag_ur(:,:,kk) = puma_ho(mask2D*angle(ur(:,:,kk)), 3,'potential',potential);%accel_unwrapping(mask2D*angle(ur(:,:,kk)),'puma',4,par.k0*simP.dn+0.025);
            else
                curr_log_imag_ur(:,:,kk) = unwrap_phase(mask2D*angle(ur(:,:,kk)));
            end
            curr_log_imag_ur(:,:,kk) = curr_log_imag_ur(:,:,kk) - mean(mean(curr_log_imag_ur(55:65,45:55,kk)));
            %imagesc(curr_log_imag_ur(:,:,kk));axis image;title(num2str(kk));colorbar;pause(0.001);
        end
        if hasandis(par,'saveUnwr')
            save(sprintf('curr_log_imag_ur_%1.2e.mat',par.rytRef),'curr_log_imag_ur');
        end
    end
end
clear roi
%%
if ~exist('curr_log_imag_ur','var')
    if exist('new_uwr','var')
        curr_log_imag_ur = new_uwr;
    elseif exist('log_imag_ur','var')
        curr_log_imag_ur = log_imag_ur;
    end
end

if size(curr_log_imag_ur,1) > size(uM.rytFBP,1)
    par.rytFBP.indx = 1 + size(curr_log_imag_ur,1)/2 - size(uM.rytFBP,1)/2:size(uM.rytFBP,1)/2 + size(curr_log_imag_ur,1)/2;
    par.rytFBP.indy = 1 + size(curr_log_imag_ur,2)/2 - size(uM.rytFBP,2)/2:size(uM.rytFBP,2)/2 + size(curr_log_imag_ur,2)/2;
    curr_log_imag_ur = curr_log_imag_ur(par.rytFBP.indx,par.rytFBP.indy,:);
else
    par.rytFBP.indx = 1 + size(uM.rytFBP,1)/2 - size(curr_log_imag_ur,1)/2:size(uM.rytFBP,1)/2 + size(curr_log_imag_ur,1)/2;
    par.rytFBP.indy = 1 + size(uM.rytFBP,2)/2 - size(curr_log_imag_ur,2)/2:size(uM.rytFBP,2)/2 + size(curr_log_imag_ur,2)/2;
    if size(Muin,3) > nnz(par.rytFBP.curr_thetas)
        uM.rytFBP = uM.rytFBP(par.rytFBP.indx,par.rytFBP.indy,par.rytFBP.curr_thetas);
        Muin = Muin(par.rytFBP.indx,par.rytFBP.indy,par.rytFBP.curr_thetas);
    end
end
if size(curr_log_imag_ur,3) > nnz(par.rytFBP.curr_thetas)
    curr_log_imag_ur = curr_log_imag_ur(:,:,par.rytFBP.curr_thetas);
end
%%
tic
par.Nryt = 0;
if exist('k_x','var')
    if hasandis(par,'Nryt')
        n_rytFBP = FBProp3D('rytov',uM.rytFBP,Muin,curr_log_imag_ur,simP,...
            par.rytFBP.thetas,k_y(par.rytFBP.curr_thetas),k_x(par.rytFBP.curr_thetas),[],par.Nryt);%simP.kz(par.rytFBP.curr_thetas));
    else
        n_rytFBP = FBProp3D('rytov',uM.rytFBP,Muin,curr_log_imag_ur,simP,...
            par.rytFBP.thetas,k_y(par.rytFBP.curr_thetas),k_x(par.rytFBP.curr_thetas),[],1);%simP.kz(par.rytFBP.curr_thetas));
    end
    %,f_rytFBP,ff_rytFBP,ff_mask]
elseif exist('kx_inc','var')
    n_rytFBP = FBProp3D('rytov',uM.rytFBP,Muin,curr_log_imag_ur,simP,...
        par.rytFBP.thetas,simP.k*sin(par.rytFBP.thetas(:,1)),...
        simP.k*sin(par.rytFBP.thetas(:,2)),[]);

else
    n_rytFBP = FBProp3D('rytov',uM.rytFBP,Muin,curr_log_imag_ur,simP,par.rytFBP.thetas);
    
end
toc
%%
if hasandis(simP,'realData') && par.negative
    n_hat.rytFBP = min(n_rytFBP,simP.n0);
else
    n_hat.rytFBP = max(n_rytFBP,simP.n0);
end

n_hat.rytFBP = setInput3D(n_hat.rytFBP,simP,simP,simP.Nroi.*ones(1,3),simP.n0);
par.rytFBP.siz = size(n_hat.rytFBP);
par.rytFBP.Lxext = simP.dx*size(n_hat.rytFBP,1); par.rytFBP.Lx_vec = simP.dx:simP.dx:par.rytFBP.Lxext;
par.rytFBP.Lyext = simP.dy*size(n_hat.rytFBP,2); par.rytFBP.Ly_vec = simP.dy:simP.dy:par.rytFBP.Lyext;
par.rytFBP.Lzext = simP.dz*size(n_hat.rytFBP,3); par.rytFBP.Lz_vec = simP.dz:simP.dz:par.rytFBP.Lzext;
par.rytFBP.dx = simP.dx; par.rytFBP.dy = simP.dy; par.rytFBP.dz = simP.dz;

uM = rmfield(uM,'rytFBP');

clear n_rytFBP;
%% n_hat.rytFBP(n_hat.rytFBP <= 1.538) = simP.n0;
myOrthoviews(n_hat.rytFBP,ceil(size(n_hat.rytFBP)/2),'Rytov',1,'w','--',1,par.negative);colormap default;
