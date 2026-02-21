%% ODT Reconstruction Born
clear
rngoi = 1;
rng(rngoi);
useGPU(0); % 1, use gpu. 0, run on cpu.
%gpuDevice(1);

% nonlinear sim measurements and born forward model 
load('./Size64/lipp_60_HalfAnglesimulation_sheppBorn_model_CG_n0_1.33_dn_0.08_Nviews_60_Nx_512_Ny_512_Nz_64_dx_0.05_dy_0.05_dz_0.05_GaussianBeam_0_photon_1_illum_circle_-42_42.mat')

% born measurements and born forward model
%load('./Size64/Born_60_HalfAnglesimulation_sheppBorn_model_CG_n0_1.33_dn_0.03_Nviews_60_Nx_512_Ny_512_Nz_64_dx_0.05_dy_0.05_dz_0.05_GaussianBeam_0_photon_1_illum_circle_-42_42.mat')
%%
uM.GTs = uM.GTc;
uems = uem;
%%
main_param_Born; % set parameters.
 
main_rytovFBP3D; % rytov approximation initial value.

input0 = (par.k*par.dz).^2*(n_hat.rytFBP.^2/simP.n0.^2 - 1);

RecoRI = @(x)sqrt((x/(par.k*par.dx)^2+1)*simP.n0^2);
input = setInput3D(input0,simP,par,par.siz,0);
BornQNLSSimBornReco;
%BornQNBornSimBornReco;