function [n,F,FF,FF_mask] = FBProp3D(approximation,uM,uin,unwrapped_imag_log_ur,simP,varargin)

lambda0 = simP.lambda0;%0.406;			% Vacuum waelength of the laser (laser diode with approx. 100 micron coherence length)
n0 = simP.n0;					% Background refractive index (immersion oil)
if isfield(simP,'NA')
    NA = simP.NA;%1.2					% Numerical aperture of the imaging objective (illumination and detection is symmetric)
else
    NA = 1.2;
end
if nargin > 5 && ~isempty(varargin{1})
    thetas = varargin{1};
    Ntheta = size(thetas,1);
else
    thetas = simP.thetas;
    Ntheta = simP.Ntheta;
end
if nargin > 6
    kx_inc = varargin{2};
    ky_inc = varargin{3};
    kz_inc = varargin{4};
end
if nargin > 9
    Niter = varargin{5};
else
    Niter = 1;
end

% Grid definition


Nx = simP.Nx;					% Number of pixels in the grid
Ny = simP.Ny;					% Number of pixels in the grid
if Nx~= size(uM,1) || Ny~= size(uM,2)
    Nx = size(uM,1);
    Ny = size(uM,2);
end
Nz = max(Nx,Ny);%simP.Nz;

fact = 1;

Nxrec = fact*Nx;
Nyrec = fact*Ny;
Nzrec = fact*Nz;
%if Nx~=Nz
%    Nz = Nx;
%end
dx = simP.dx;%camera_pixel_size/magnification;	% Effective pixel size in the object plane
dy = simP.dy;%camera_pixel_size/magnification;	% Effective pixel size in the object plane
dz = simP.dz;

Lx = dx*Nx;					% Width of the computation window
Ly = dy*Ny;					% Width of the computation window
Lz = dz*Nz;

dkx = 2*pi/Lx;				% Pixel size in the Fourier domain
dky = 2*pi/Ly;				% Pixel size in the Fourier domain

dkxrec = 2*pi/(Nxrec*dx);
dkyrec = 2*pi/(Nyrec*dy);
dkzrec = 2*pi/(Nzrec*dz);

x = (-Nx/2:Nx/2-1)*dx;		% Coordinate vector in the spatial domain
y = (-Ny/2:Ny/2-1)*dy;		% Coordinate vector in the spatial domain
kx = (-Nx/2:Nx/2-1)*dkx;	% Coordinate vector in the Fourier domain
ky = (-Ny/2:Ny/2-1)*dky;	% Coordinate vector in the Fourier domain



[YY,XX] = meshgrid(y,x);

[Kyy,Kxx] = meshgrid(ky,kx);

K = sqrt(Kxx.^2 + Kyy.^2);

k0 = 2*pi/lambda0;			% Vaccum k-vector
k = k0*n0;					% k-vector in the background medium of refractive index n0
kmax = k0*NA;				% Largest k-vector transmitted by an optical system of given NA
objective_aperture_ind = find(K < kmax);	% Array containing the indices of the pixels within the aperture of radius kmax
Nind = length(objective_aperture_ind);				% Number of pixels within the aperture of radius kmax

% Incident plane waves vectors

if ~exist('kx_inc','var')
    kx_inc = k*sin(thetas(:,2));
    ky_inc = k*sin(thetas(:,1));
end

if ~exist('kz_inc','var') || isempty(kz_inc)
    kz_inc = sqrt(k.^2 - kx_inc.^2 - ky_inc.^2);	% Because the object is immersed in background index n0, we have kx^2 + ky^2 + kz^2 = k^2
end

if iscolumn(kx_inc)
    s0 = [kx_inc ky_inc kz_inc]/k;		% Normalized k-vectors of the incident waves
else
    s0 = [kx_inc' ky_inc' kz_inc']/k;		% Normalized k-vectors of the incident waves
end
sx = Kxx(objective_aperture_ind)/k;		% Normalized k-vectors of the detection directions. The minus sign is because of the orientation of the detection plane
sy = Kyy(objective_aperture_ind)/k;
sz = real(sqrt(1 - sx.^2 - sy.^2));
s = [sx sy sz];


FF = zeros(Nxrec, Nyrec, Nzrec);							% FF will contain the 3D Fourier transform of the reconstructed scattering potential F
FF_mask = false(Nxrec, Nyrec, Nzrec);					% Set of the voxels that will be populated in F (for monitoring only)
cumulated_contributions = zeros(Nxrec, Nyrec, Nzrec);	% Arrays in which we put the counts in case a voxel gets assigned more than once


for kk = 1:Ntheta
    u = uM(:,:,kk);
    ui = uin(:,:,kk);
    ur = u./ui;
    us = u - ui;
    
    if strcmpi(approximation, 'iWolf')
        um_0 = ur;%should be amplitude only
    elseif strcmpi(approximation, 'Born')
        % Born approximation
        um_0 = us./ui;
    elseif strcmpi(approximation, 'Rytov')
        % Rytov approximation
        imag_log_ur = unwrapped_imag_log_ur(:,:,kk) - mean(mean(unwrapped_imag_log_ur(:,:,kk)));%imag(log(ur));
        real_log_ur = real(log(ur));
        
        imag_log_ur(isnan(imag_log_ur)) = 0;
        real_log_ur(isnan(real_log_ur)) = 0;
        um_0 = real_log_ur + 1i*imag_log_ur;
    else
        error('Invalid approximation')
    end
    
    um_0 = um_0 .* exp(1i*(kx_inc(kk)*XX + ky_inc(kk)*YY));	% We modulate the field by the incident plane wave. Note that the field are initially corrected for the illumination plane wave so that we have to put it back here
    far_field = fftshift(fft2(fftshift(um_0)));
    far_field = k*sz/(1i*2*pi*dz).*far_field(objective_aperture_ind);
    
    S = s - repmat(s0(kk, :), Nind, 1);		% According to Wolf formula, S is the Fourier component that is probed by wave s under illumination s0, with S = s - s0
    K_nearest_neighbor = round(S./repmat([dkxrec dkyrec dkzrec]/k, Nind, 1));	% We calculate the normalized S vector in pixel units. That will be the surface of the Ewald sphere. We do simple nearest neighbor interpolation
    px = K_nearest_neighbor(:, 1) + round((Nxrec-1)/2)+1;                  % More compact notation and set the center of the coordinate system to the corner of the array
    py = K_nearest_neighbor(:, 2) + round((Nyrec-1)/2)+1;
    pz = K_nearest_neighbor(:, 3) + round((Nzrec-1)/2)+1;
    
    fourier_space_indices = (pz-1)*Nxrec*Nyrec + (py-1)*Nxrec + px;	% Calculate the indices in the Fourier domain
    
    FF(fourier_space_indices) = FF(fourier_space_indices) + far_field;	% Assign the far field data to the surface of the Ewald sphere
    
    FF_mask(fourier_space_indices) = true;	% Set FF_mask to 1 on the Ewald sphere
    
    cumulated_contributions(fourier_space_indices) = cumulated_contributions(fourier_space_indices) + 1;	% We count an increment of 1 for each voxel that we have just populated
    
    
end


FF = fact^2*FF./max(cumulated_contributions, 1);	% We divide the populated Fourier spectrum by the cumulated contributions

doIter = Niter > 1;
if doIter
    %max_b = 0.4/(4*pi)*k0^2;
    FFmeas = FF(FF_mask);
    figure(99);
    for jj = 1:Niter
        F = fftshift(ifftn(ifftshift(FF)));	% We reconstruct the scattering potential defined as F = k^2/(4pi)(n^2 - n0^2)
        F = max(real(F),0) + 1i*imag(F);%Positivity constraint
        %F = min(real(F),max_b) + 1i*imag(F);
        %F(real(F) > max_b) = 1i*imag(F(real(F) > max_b));
        FF = fftshift(fftn(ifftshift(F)));
        FF(FF_mask) = FFmeas;
        subplot(121);imagesc(log(1 + abs(FF(:,:,round(end/2)))));axis image;
        subplot(122);imagesc(real(sqrt(4*pi*F(:,:,round(end/2))/k0^2 + n0^2)));colorbar;axis image;
        %suptitle(sprintf('Iteration %i',jj));drawnow
    end
end
F = fftshift(ifftn(ifftshift(FF)));	% We reconstruct the scattering potential defined as F = k^2/(4pi)(n^2 - n0^2)
if doIter
    F = max(real(F),0) + 1i*imag(F);%Positivity constraint
end
%F = min(real(F),max_b) + 1i*imag(F);
%F(real(F) > max_b) = 1i*imag(F(real(F) > max_b));
n = real(sqrt(4*pi*F/k0^2 + n0^2));	% Retrieve the refractive index from the scattering potential

end


