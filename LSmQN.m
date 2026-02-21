%% MAIN 3D FastLipp

fpath = createFolder(par.folder,'LSSimLSReco',par.param_fname);

fprintf('%1.2f < 0.5 ? \n', par.dx/simP.lambda0*simP.n0);
%% Set 3D incident fields

cuin = curr_uem;

uin_lipp = setUin3D(par,simP,'lipp',single(cuin),SP);

%% Set Measurements
if par.discreteGtild
    if par.SP.nCamera > 1
        uem_lipp = curr_uem(:,:,:,par.curr_thetas);
    else
        uem_lipp = curr_uem(:,:,par.curr_thetas);
    end
else
    if exist('curr_uem','var')
        uem_lipp = curr_uem(:,:,par.curr_thetas);
    elseif ~isempty(cuin)
        uem_lipp = cuin(:,:,par.curr_thetas);
    else
        error('Please provide Uem or Uin at the sensor positions');
    end
end

[uM.lipp,M_lipp,uem_lipp] = setMeasurements3D(uM.curr_y,par,simP,par.noreflection,'lipp',...
    M_curr,curr_SP, uem_lipp);

%% Set the forward model and regularization

H = setmultOp3Dnew(uin_lipp,M_lipp,par,'lipp',uem_lipp,curr_SP, simP);

F = createCostsArray(H,uM.lipp,true);
if 0
    F = cell(par.Ntheta,1);
    repel = repmat({':'},1,max(ndims(uM.lipp)-1,1));
    for kk = 1:par.Ntheta
        F{kk} = CostL1(H.mapsCell{kk}.sizeout, ...
            uM.lipp(repel{:},kk)/squeeze(sum(sum(abs(uM.lipp(repel{:},kk)).^2,1),2)))...
            *H.mapsCell{kk};
    end
    F = CostSummation(F,1);
end

F = F.makePartialSummation(par.Lsub);
F.setPartialGrad(par.stochGrad);

switch par.RegType
    case 'TV'
        if isfield(par,'mu') && par.mu <= 1
            Reg = CostMixNorm21([par.siz,3],4)*LinOpGrad(par.siz,[],'zeros',1./[par.mu,par.mu,1 - par.mu]);
        else
            Reg = CostMixNorm21([par.siz,3],4)*LinOpGrad(par.siz,[],'zeros');
        end
    case 'Hess'
        Hess = LinOpHess(par.siz,'circular');
        Reg = CostMixNormSchatt1(Hess.sizeout,par.p)*Hess;
        Reg.setC(par.bounds(1),par.bounds(2));
        Reg.max_iter = par.prox_iter;
        Reg.xtol = par.xtolin;
end

%%
if exist('n','var') && ~isempty(n)
    n_gt.lipp = setInput3D(n,simP,par,par.siz,simP.n0);
    f_res.lipp = (par.k*par.dz)^2*(n_gt.lipp.^2/simP.n0^2 - 1);

    figure(10);dispGT3D(par,simP,[],f_res.lipp,par.negative,[],gcf);
else
    f_res.lipp = [];
    n_gt.lipp = zeros(par.siz);
end
if strcmp(par.RegType,'TV')
    Reg.setProxAlgo(par.bounds,par.prox_iter,par.xtolin,...
        myOutputOpti(1,n_gt.lipp,5,par.k,simP.n0,'Lipp',par.negative));
    Reg.optim.ItUpOut = 0;
end

%% Test forward
aroi = min(10,par.Ntheta);
tic
if isempty(f_res.lipp)
    u_gt = gather(H.mapsCell{aroi}*(gpuCpuConverter(input)));
else
    u_gt = gather(H.mapsCell{aroi}*(gpuCpuConverter(f_res.lipp)));
end
%compareFields(u_gt, uM.GTc(:,:,aroi) - uem(:,:,aroi));
compareFields(u_gt, uM.lipp(:,:,aroi));
drawnow;
toc

%====================================
% compute the cost of the initial value
temp_cost = zeros(1,3);
temp_aa = F;
Reg_temp = Reg;
temp_aa.setPartialGrad(0)
temp_cost(1,1) = gather(temp_aa.apply(gpuCpuConverter(input))); % check how to compute the cost of the data fidelty
temp_cost(1,2) = gather(Reg_temp.apply(gpuCpuConverter(input)));
%temp_cost(1,3) = gather(Reg_temp.apply(gpuCpuConverter(input)));
%====================================
%% Run FISTA for Lippmann
save(fullfile(fpath,'F.mat'),'F','-v7.3');
save(fullfile(fpath,'RecoOperator.mat'),'RecoRI');
%%
% test the choice of parameters -- the number of subset & and choice of

par.saveX = 1;      % 1; 100; %3;
par.ItUpOut = 1;    % 1; 3;
MaxIter = 120;
par.prox_iter = 200;


Numsubset = 15;
algName_set = {'QN_ours_SR1_UR'};%'FBP', 'QN_ours_SR1','QN_ours_full_SR1','FBP_full'

isComputed = true;

for ll = 1:length(algName_set)
    algName = algName_set{ll};

    for kk = 1:length(regul_set)
        par.regul = regul_set(kk);
        for JJ = 1:length(Numsubset)
            par.Lsub = Numsubset(JJ);
            OutOp = myOutputOpti(1,f_res.lipp,par.iterVerb,[par.k*par.dz,simP.k],simP.n0,'lipp',ceil(par.siz/2));
            OutOp.doPlot = false;
            OutOp.saveX = par.saveX;
            if isfield(par,'saveAll')
                OutOp.saveAll = par.saveAll;
            end

            F.counterSet = false;
            F.setPartialGrad(2);
            F.setLsub(par.Lsub);
            if isa(Reg,'CostRectangle')
                %F.alpha(end) = par.regul;
            elseif isa(Reg,'CostMultiplication')
                Reg = par.regul*Reg.cost2;
            else
                Reg = par.regul*Reg;
            end

            switch algName

                case 'FBP'
                    opt = OptiFBS(F, Reg);
                    opt.fista = strcmp(par.accel,'fista');
                    opt.gam = par.gam; %100 Linear, has a Lipschitz
                    opt.OutOp = OutOp;

                case 'QN_ours_SR1'

                    if isa(Reg,'CostMultiplication')
                        H2 = Reg.cost2.H2;
                    else
                        H2 = Reg.H2;
                    end

                    opt = OptiBatchQNP(F,H2,par.regul,'iso',par.bounds,OutOp);
                    opt.UpdateHessianType = 'SRI';
                    opt.isComputerErr = isComputed; % true, compute the gradient error.
                    opt.gamma = 0.8;
                    opt.a_k = 1;
                    opt.maxiterin = par.prox_iter;
                    opt.L = 1/par.gam;
                case 'QN_ours_SR1_UR'

                    if isa(Reg,'CostMultiplication')
                        H2 = Reg.cost2.H2;
                    else
                        H2 = Reg.H2;
                    end

                    opt = OptiBatchQNPUR(F,H2,par.regul,'iso',par.bounds,OutOp);
                    
                    opt.UpdateHessianType = 'SRI';
          
                    opt.isAnchor = true;
                    opt.isAnchorGrad = true;
                    opt.isComputerErr = isComputed;
                    opt.gamma = 0.8;
                    opt.a_k = 1;
                    opt.maxiterin = par.prox_iter;
                    opt.L = 1/par.gam;
                case 'QN_ours_full_SR1'
                    par.Lsub = par.Ntheta;
                    F.setLsub(par.Lsub);
                    if isa(Reg,'CostMultiplication')
                        H2 = Reg.cost2.H2;
                    else
                        H2 = Reg.H2;
                    end
                    opt = OptiBatchQNP(F,H2,par.regul,'iso',par.bounds,OutOp);
                    opt.UpdateHessianType = 'SRI';
                    opt.gamma = 0.8;
                    opt.a_k = 1;
                    opt.maxiterin = par.prox_iter;
                    opt.L = 1/par.gam;

                case 'FBP_full'
                    par.Lsub = par.Ntheta;
                    F.setLsub(par.Lsub);
                    opt = OptiFBS(F, Reg);
                    opt.fista = strcmp(par.accel,'fista');
                    opt.gam = par.gam;
                    opt.OutOp = OutOp;
                otherwise
                    error('Alg. is not defined.\n')
            end

            opt.ItUpOut = par.ItUpOut; % save the results at each ItUpOut iteration.
            opt.CvOp = TestCvgStepRelative(par.xtol);
            opt.maxiter = MaxIter;
            opt.run(gpuCpuConverter(input));

            n_hat.lipp = gather(padarray(sqrt(opt.xopt/(par.k*par.dz)^2 + 1)*simP.n0,...
                (par.siz - par.siz)/2, simP.n0));
            fname = sprintf('%s',simP.ObjectType);
            if strcmp(simP.ObjectType,'mult_bead') && isfield(simP,'nBeads')
                fname = sprintf('%s_%i_low',fname,simP.nBeads);
            end
            if isfield(simP,'dn')
                fname = sprintf('%s_dn_%1.2f',fname,simP.dn(end));
            end
            if isempty(n_gt.lipp)
                curr_snr = nan;
            else
                curr_snr = snr(n_gt.lipp,n_gt.lipp - n_hat.lipp);
            end

            switch algName

                case {'QN_ours_SR1_UR','QN_ours_SR1','QN_ours_full_SR1'}
                    fname = sprintf('%sGamma_%1.2f_Step_%1.2f_Subsetno_%i_%s_lipp_regul_%1.2e_Nx_%i_Nxext_%i_Nz_%i_Nzext_%i_iterin_%i_snr_%1.2f_gam_%1.2e_RegType_%s.mat',...
                        algName,opt.gamma,opt.a_k,par.Lsub,fname, par.regul,par.Nx,par.Nxext,par.Nz,par.Nzext,...
                        opt.maxiterin,curr_snr,par.gam,par.RegType);
                case {'FBP','FBP_full'}
                    fname = sprintf('%sStep_%1.2f_Subsetno_%i_%s_lipp_regul_%1.2e_Nx_%i_Nxext_%i_Nz_%i_Nzext_%i_iterin_%i_snr_%1.2f_gam_%1.2e_RegType_%s.mat',...
                        algName,par.gam,par.Lsub,fname,par.regul,par.Nx,par.Nxext,par.Nz,par.Nzext,...
                        par.Niterin,curr_snr,par.gam,par.RegType);
            end

            fprintf('Saving data as %s...\n',fname);
            OutOp = opt.OutOp;
            curr_nhat = n_hat;
            curr_uM = gather(uM.lipp);
            curr_snr = gather(curr_snr);
            n_hat = struct('lipp',curr_nhat.lipp);
            if isComputed
                grad_err_set = gather(opt.grad_err_set);
                yk_err_set = gather(opt.yk_err_set);
                grad_norm_full = gather(opt.grad_norm_full);
                grad_norm_estimate = gather(opt.grad_norm_estimate);
                save(fullfile(fpath,strcat(sprintf('Comp%i_',isComputed),fname)),'OutOp','grad_err_set','yk_err_set','grad_norm_full','grad_norm_estimate','input','par','n_hat','curr_uM','Reg','temp_cost','-v7.3');
            else
                save(fullfile(fpath,strcat(sprintf('Comp%i_',isComputed),fname)),'OutOp','input','par','n_hat','curr_uM','Reg','temp_cost','-v7.3');
            end
            %save(fullfile(fpath,strcat(sprintf('Comp%i_',isComputed),fname)),'OutOp','input','par','n_hat','curr_uM','Reg','temp_cost','-v7.3');
            %n_hat  curr_nhat;
            clear curr_nhat
        end
    end
end
%%
%{
% check the choice of gamma
gamma_set = [0.6;0.7;0.8;0.85;0.9;1];
par.Lsub = 15;
%algName_set = {'FBP','QN_ours_SR1','QN_anchorIter_SR1','QN_anchorGrad_SR1','QN_anchorIterGrad_SR1'};
algName_set = {'QN_ours_SR1_UR','QN_anchorIterGrad_SR1'};
isComputed = false;
for ll = 1:length(algName_set)
    algName = algName_set{ll};
    for kk = 1:length(regul_set)%[1:2:10,2:2:10]
        %Reg = Reg_temp;
        par.regul = regul_set(kk);

        for JJ =1:length(gamma_set)

            OutOp = myOutputOpti(1,f_res.lipp,par.iterVerb,[par.k*par.dz,simP.k],simP.n0,'lipp',ceil(par.siz/2));%[135,112,50]);%
            OutOp.doPlot = false; %par.doPlot;
            OutOp.saveX = par.saveX;
            if isfield(par,'saveAll')
                OutOp.saveAll = par.saveAll;
            end

            F.counterSet = false;
            F.setPartialGrad(2);
            F.setLsub(par.Lsub);

            if isa(Reg,'CostRectangle')
                %F.alpha(end) = par.regul;
            elseif isa(Reg,'CostMultiplication')
                Reg = par.regul*Reg.cost2;
            else
                Reg = par.regul*Reg;
            end


            if isa(Reg,'CostMultiplication')
                H2 = Reg.cost2.H2;
            else
                H2 = Reg.H2;
            end

            if gamma_set(JJ)==0.8
                continue;
            end

            switch algName

                case {'QN_ours_SR1Computed','QN_ours_SR1','QN_anchorIter_SR1','QN_anchorGrad_SR1','QN_anchorIterGrad_SR1','QN_ours_ModifiedSR1','QN_anchorIter_ModifiedSR1','QN_anchorGrad_ModifiedSR1','QN_anchorIterGrad_ModifiedSR1'}

                    if isa(Reg,'CostMultiplication')
                        H2 = Reg.cost2.H2;
                    else
                        H2 = Reg.H2;
                    end

                    opt = OptiWeightedProx_Tao20241230(F,H2,par.regul,'iso',par.bounds,OutOp);
                    if strcmp(algName,'QN_ours_ModifiedSR1') || strcmp(algName,'QN_anchorIter_ModifiedSR1') || strcmp(algName,'QN_anchorGrad_ModifiedSR1') || strcmp(algName,'QN_anchorIterGrad_ModifiedSR1')
                        opt.UpdateHessianType = 'ModifiedSRI';
                    else
                        opt.UpdateHessianType = 'SRI';
                    end

                    if strcmp(algName,'QN_anchorIter_ModifiedSR1') || strcmp(algName,'QN_anchorIterGrad_ModifiedSR1') || strcmp(algName,'QN_anchorIter_SR1') || strcmp(algName,'QN_anchorIterGrad_SR1')
                        opt.isAnchor = true;
                    else
                        opt.isAnchor = false;
                    end

                    if strcmp(algName,'QN_anchorGrad_ModifiedSR1') || strcmp(algName,'QN_anchorIterGrad_ModifiedSR1') || strcmp(algName,'QN_anchorGrad_SR1')  ||  strcmp(algName,'QN_anchorIterGrad_SR1')
                        opt.isAnchorGrad = true;
                    else
                        opt.isAnchorGrad = false;
                    end

                    opt.isComputerErr = isComputed; % true, compute the gradient error.
                    opt.gamma = gamma_set(JJ);
                    opt. ALGType = 'Type-II';
                    opt.a_k = 1;
                    opt.maxiterin = par.prox_iter;
                    opt.L = 1/par.gam;
                case {'QN_ours_SR1_UR','QN_ours_ModifiedSR1_UR'}

                    if isa(Reg,'CostMultiplication')
                        H2 = Reg.cost2.H2;
                    else
                        H2 = Reg.H2;
                    end

                    opt = OptiWeightedProx_Tao20250113UR(F,H2,par.regul,'iso',par.bounds,OutOp);
                    if strcmp(algName,'QN_ours_SR1_UR')
                        opt.UpdateHessianType = 'SRI';
                    elseif strcmp(algName,'QN_ours_ModifiedSR1_UR')
                        opt.UpdateHessianType = 'ModifiedSRI';
                    end

                    opt.isAnchor = true;
                    opt.isAnchorGrad = true;

                    opt.isComputerErr = isComputed;
                    opt.gamma = gamma_set(JJ);
                    opt. ALGType = 'Type-II';
                    opt.a_k = 1;
                    opt.maxiterin = par.prox_iter;
                    opt.L = 1/par.gam;

                case  {'QN_ours_full_SR1','QN_ours_full_ModifiedSR1'}
                    par.Lsub = par.Ntheta;
                    F.setLsub(par.Lsub);
                    if isa(Reg,'CostMultiplication')
                        H2 = Reg.cost2.H2;
                    else
                        H2 = Reg.H2;
                    end
                    opt = OptiWeightedProx_Tao_Mar2022(F,H2,par.regul,'iso',par.bounds,OutOp);
                    if strcmp(algName,'QN_ours_full_SR1')
                        opt.UpdateHessianType = 'SRI';
                    else
                        opt.UpdateHessianType = 'ModifiedSRI';
                    end
                    opt.gamma = gamma_set(JJ);
                    opt. ALGType = 'Type-II';
                    opt.a_k = 1;
                    opt.maxiterin = par.prox_iter;
                    opt.L = 1/par.gam;

                otherwise
                    error('Alg. is not defined.\n')
            end

            opt.ItUpOut = par.ItUpOut; % save the results at each ItUpOut iteration.
            opt.CvOp = TestCvgStepRelative(par.xtol);
            opt.maxiter = MaxIter;
            opt.run(gpuCpuConverter(input));

            n_hat.lipp = gather(padarray(sqrt(opt.xopt/(par.k*par.dz)^2 + 1)*simP.n0,...
                (par.siz - par.siz)/2, simP.n0));
            fname = sprintf('%s',simP.ObjectType);
            if strcmp(simP.ObjectType,'mult_bead') && isfield(simP,'nBeads')
                fname = sprintf('%s_%i_low',fname,simP.nBeads);
            end
            if isfield(simP,'dn')
                fname = sprintf('%s_dn_%1.2f',fname,simP.dn(end));
            end
            if isempty(n_gt.lipp)
                curr_snr = nan;
            else
                curr_snr = snr(n_gt.lipp,n_gt.lipp - n_hat.lipp);
            end

            switch algName

                case {'QN_ours_SR1_UR','QN_ours_ModifiedSR1_UR','QN_ours_SR1Computed','QN_ours_SR1','QN_anchorIter_SR1','QN_anchorGrad_SR1','QN_anchorIterGrad_SR1','QN_ours_ModifiedSR1','QN_anchorIter_ModifiedSR1','QN_anchorGrad_ModifiedSR1','QN_anchorIterGrad_ModifiedSR1','QN_wang','QN_ours_full_SR1','QN_ours_full_ModifiedSR1'}
                    fname = sprintf('%sGamma_%1.2f_Step_%1.2f_Subsetno_%i_%s_lipp_regul_%1.2e_Nx_%i_Nxext_%i_Nz_%i_Nzext_%i_iterin_%i_snr_%1.2f_gam_%1.2e_RegType_%s.mat',...
                        algName,opt.gamma,opt.a_k,par.Lsub,fname, par.regul,par.Nx,par.Nxext,par.Nz,par.Nzext,...
                        opt.maxiterin,curr_snr,par.gam,par.RegType);
                case {'FBP','FBP_full'}
                    fname = sprintf('%sStep_%1.2f_Subsetno_%i_%s_lipp_regul_%1.2e_Nx_%i_Nxext_%i_Nz_%i_Nzext_%i_iterin_%i_snr_%1.2f_gam_%1.2e_RegType_%s.mat',...
                        algName,par.gam,par.Lsub,fname,par.regul,par.Nx,par.Nxext,par.Nz,par.Nzext,...
                        par.Niterin,curr_snr,par.gam,par.RegType);
            end

            fprintf('Saving data as %s...\n',fname);
            OutOp = opt.OutOp;
            curr_nhat = n_hat;
            curr_uM = gather(uM.lipp);
            curr_snr = gather(curr_snr);
            n_hat = struct('lipp',curr_nhat.lipp);
            save(fullfile(fpath,strcat(sprintf('Comp%i_',isComputed),fname)),'OutOp','input','par','n_hat','curr_uM','Reg','temp_cost','-v7.3');
            clear curr_nhat
        end
    end
end
%}