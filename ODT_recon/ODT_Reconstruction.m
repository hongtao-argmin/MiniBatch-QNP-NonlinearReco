%% ODT Reconstruction
close all;
clear all;
%reset(gpuDevice(2))
rng(1)
%gpuDevice(2);
useGPU(0);

% load the measurements.
foldername = './Size64/';
filename_set = {
'lipp_60_HalfAnglesimulation_shepp_model_CG_n0_1.33_dn_0.09_Nviews_60_Nx_512_Ny_512_Nz_64_dx_0.05_dy_0.05_dz_0.05_GaussianBeam_0_photon_1_illum_circle_-42_42.mat';
  };

%%
for itername_index = 1:length(filename_set)
    close all;
    %reset(gpuDevice(2))
    rng(1)
    %gpuDevice(2);
    useGPU(0);
    load([foldername filename_set{itername_index}])
simP.ObjectType = sprintf('%s_%1.2f',simP.ObjectType,simP.dn(1));
%% Load simulated data
fac = 1;
if ~isfield(simP,'NA')
	simP.NA = simP.n0;
end
uM.GTs = uM.GTc;
uems = uem;
%%
main_param_LSm; % Set the parameters
par.simP = simP;
main_rytovFBP3D;	
input0 = (par.k*par.dx)^2*(n_hat.rytFBP.^2/simP.n0^2 - 1);
RecoRI = @(x)sqrt((x/(par.k*par.dx)^2+1)*simP.n0^2);
input = setInput3D(input0,simP,par,par.siz,0);
LSmQN;
end