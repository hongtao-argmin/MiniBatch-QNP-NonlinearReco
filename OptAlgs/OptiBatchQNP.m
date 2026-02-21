classdef OptiBatchQNP < Opti
    % OptiWeightedProx_Tao2024NewTest
    % **************************
    % Tao:
    % Implement Type-I and II QNP for problems with data
    % fidelity + TV while the data fidelity term has multiple measurements.
    % Test the error of several criteria and diff ways to estimate the Hessian matrix
    % **************************
    % Protected Set and public Read properties
    properties (SetAccess = protected,GetAccess = public)
        P;
        reg;
        D;
        ndims;
        bounds;
        CostR;
        sizein;

        % for defining block version.
        x_k_cell;
        g_k_cell;
        y_k_cell
        Hessian_u_k_cell;
        Hessian_D_k_cell;
        Numblock;
        grad_err_set;
        yk_err_set;
        grad_norm_full;
        grad_norm_estimate;
        F_temp;
        options = optimoptions('fsolve','Display','off');%,'iter','Tolx');
    end
    % Full public properties
    properties

        a_k = 0.95;            % stepsize for the quasi-newton / we may do line search if this does not work.
        gamma = 0.8;           % overshot parameter, should be small than 1.
        UpdateHessianType = 'SRI';
        isAnchor = false;      % true, we always use x_k rather than the aggregated whole previous histories
        isAnchorGrad = false;  % true, we only use \nabla f(x_k) rather than the aggregated whole previous histories
        isComputerErr = false; % true, compute the gradient error.
        L;
        xtolTV = 1e-6;
        maxiterin = 200;
        TVtype = 'iso'; %'l1' or 'iso'
        F;
        MaxIterProj;
    end

    methods
        %% Constructor
        function this = OptiBatchQNP(F,D,reg,TVtype,bounds,OutOp)
            if nargin < 6 || isempty(OutOp), OutOp = OutputOpti(1);end
            if nargin < 5 || isempty(bounds), bounds = [-inf,inf];end
            this.name = 'OptiQNP';
            this.F = F;
            this.F_temp = F;
            this.D = D;
            this.ndims = length(D.sizein);
            this.reg = reg;
            this.TVtype = TVtype;
            this.bounds = bounds;
            this.CostR = CostRectangle(D.sizein,bounds(1),bounds(2));
            if strcmp(this.TVtype,'iso')
                this.cost= F + reg*CostMixNorm21(D.sizeout,length(D.sizeout))*this.D;
            elseif strcmp(this.TVtype,'l1')
                this.cost= F + reg*CostL1(D.sizeout,length(D.sizeout))*this.D;
            else
                error('Unknown type of TV');
            end
            this.OutOp=OutOp;
            this.L = F.lip;
        end

        function initialize(this,x0)
            % Reimplementation from :class:`Opti`.
            assert(this.L > 0,'L must be larger than 0');
            this.sizein = size(x0);
            this.MaxIterProj = 100;
            this.Numblock = ceil(this.F.numMaps/this.F.Lsub);
            this.x_k_cell = cell(this.Numblock,1);
            this.g_k_cell = cell(this.Numblock,1);
            this.y_k_cell = cell(this.Numblock,1);
            this.grad_err_set = zeros_(this.maxiter,1);
            this.yk_err_set = zeros_(this.maxiter,1);
            this.grad_norm_full = zeros_(this.maxiter,1);
            this.grad_norm_estimate = zeros_(this.maxiter,1);
            this.Hessian_u_k_cell = cell(this.Numblock,1);
            this.Hessian_D_k_cell = zeros_(this.Numblock,1);
            initialize@Opti(this,x0);
        end

        function flag = doIteration(this)
            % evaluate full gradient
            if this.isComputerErr
                if this.niter>=this.Numblock
                    g_k_total = zeros_(this.sizein);
                    for iter_index = 1:this.Numblock
                        temp_grad = this.F_temp.applyGrad(this.xopt);
                        g_k_total = g_k_total+temp_grad;
                    end
                    g_k_total =  g_k_total/this.Numblock;
                end
            end
            % evaluate gradient
            g_k_1 = this.F.applyGrad(this.xopt);
            this.xold = this.xopt;
            indoi = mod(this.niter,this.Numblock)+1;
            % first iteration // only need to formulate B^{-1}
            if  this.niter == 0 || this.niter< this.Numblock%mod(this.niter,30) == 0 || (mod(this.niter,30) < this.Numblock)%
                % run ISTA for the beginning iterations.
                this.x_k_cell{indoi} = this.xopt;
                this.g_k_cell{indoi} = g_k_1;
                D_inv = 1/(this.L*this.a_k)*ones_(1);
                D_diag = 1/D_inv;
                u_inv = zeros_(this.sizein);
                u = zeros_(this.sizein);
            else
                % estimate Hessian
                s_k = this.xopt - this.x_k_cell{indoi};
                y_k = g_k_1 - this.g_k_cell{indoi};

                this.x_k_cell{indoi} = this.xopt;
                this.g_k_cell{indoi} = g_k_1;
                if this.isComputerErr
                    this.y_k_cell{indoi} = y_k;
                    count_y_k_temp = 0;
                    % compute the gradient error full vs approximated
                    for iter=1:this.Numblock
                        if iter==1
                            temp_grad_est = this.g_k_cell{iter};
                            if isempty(this.y_k_cell{iter})
                                temp_y_k = 0;
                            else
                                temp_y_k = this.y_k_cell{iter};
                                count_y_k_temp = count_y_k_temp+1;
                            end
                        else
                            temp_grad_est = temp_grad_est+this.g_k_cell{iter};
                            if ~isempty(this.y_k_cell{iter})
                                temp_y_k = temp_y_k+this.y_k_cell{iter};
                                count_y_k_temp = count_y_k_temp+1;
                            end
                        end
                    end
                    temp_grad_est = temp_grad_est/this.Numblock;
                    this.grad_err_set(this.niter+1) = norm(temp_grad_est-g_k_total,'fro');
                    this.grad_norm_full(this.niter+1) = norm(g_k_total,'fro');
                    this.grad_norm_estimate(this.niter+1) = norm(temp_grad_est,'fro');
                    this.yk_err_set(this.niter+1) = norm(temp_y_k/count_y_k_temp,'fro');
                end
                % H = D_diag+ sign_u uu^T, sign_u is set to be negative
                [D_diag,u,D_inv,u_inv] = this.EstimateHessian(s_k,y_k);
            end

            this.Hessian_u_k_cell{indoi} = u;      % store vector u.
            this.Hessian_D_k_cell(indoi) = D_diag; % store the scale variable which is diagonal.
            
            if this.niter<this.Numblock
                Y = this.xopt-1/this.L*g_k_1;
            else
                % we begin to accumulate the histories after first round.
                % implement \sum_(B_i) B_i * x_i
                Y = this.Hessian_D_k_cell(1)*this.xopt-...
                            this.Hessian_u_k_cell{1}*(sum(this.Hessian_u_k_cell{1}.*this.xopt,'all'))-this.a_k*g_k_1;

                for index_Hessian = 2:this.Numblock
                    if isempty(this.Hessian_u_k_cell{index_Hessian})
                        continue;
                    elseif any(this.Hessian_u_k_cell{index_Hessian}(:))
                        Y = Y+this.Hessian_D_k_cell(index_Hessian)*this.xopt-...
                                    this.Hessian_u_k_cell{index_Hessian}*...
                                    (sum(this.Hessian_u_k_cell{index_Hessian}.*this.xopt,'all'))-...
                                    this.a_k*g_k_1;
                    else
                        Y = Y+this.Hessian_D_k_cell(index_Hessian)*this.xopt-...
                                    this.a_k*g_k_1;
                    end
                end
            end

            if this.niter<this.Numblock
                this.xopt = this.ProxTVWarm(Y,D_inv,u_inv);
            else
                this.xopt = this.ProxTVWarm_NV(Y);
            end
            flag = this.OPTI_NEXT_IT;
        end

        function [D_diag,u,D_inv,u_inv] = EstimateHessian(this,s_k,y_k)
            % s_k: diff x_k
            % y_k: diff \nabla f(x+k)
            % estimate Hessian
            % D_inv,u_inv,sign_u_inv
            % the whole function estimates the inverse of Hessian matrix but the return will be the Hessian matrix
            % one can change to a better Hessian estimation can improve the
            % performance. 
            % e.g., ref. -> Tao Hong et al., Convergent Complex Quasi-Newton Proximal Methods for Gradient-Driven Denoisers in Compressed Sensing MRI Reconstruction, IEEE Transactions on Computational Imaging, vol. 11, pp. 1534-1547, Oct. 2025
            % that can guarantee the positive definite of the Hessian even
            % the function is nonconvex. We also do not need the bounded
            % Hessian assumption presented in the paper.
            switch this.UpdateHessianType
                case 'SRI'
                    % estimate from inverse H^{-1}
                    % one can also estimate the Hessian from the H directly
                    % 
                    tau_BB2 = sum(s_k.*y_k,'all')/sum(y_k.^2,'all');
                    if tau_BB2<0
                        D_inv = this.a_k/(this.L)*ones_(1);
                        D_diag = 1/D_inv;
                        u_inv = zeros_(this.sizein);
                        u = zeros_(this.sizein);
                    else
                        H_0 = tau_BB2*this.gamma;
                        temp_1 = s_k-H_0*y_k;
                        temp_2 = sum(temp_1.*y_k,'all');
                        D_inv = H_0;
                        D_diag = 1/D_inv;
                        if temp_2<=1e-8*norm(y_k(:))*norm(temp_1(:))
                            u_inv = zeros_(this.sizein);
                            u = zeros_(this.sizein);
                        else
                            u_inv = temp_1/sqrt(temp_2);
                            u = D_inv\u_inv/sqrt(1+sum(u_inv.*(D_inv\u_inv),'all'));
                        end
                    end
                otherwise
                    % put more strategies of estimating Hessian or weighting
                    % here. or design a fixed preconditioner.
                    error([this.UpdateHessianType ' is not defined.\n'])
            end
        end

        function v_out = ProxTVWarm_NV(this,v)
            % Hessian matrix B is B = D \pm u*u'; -- test on positive first.
            % -------------------------------------------------------------------------
            alph = this.a_k*this.reg*this.Numblock;
            V_acc = zeros_(prod(this.sizein),this.Numblock);
            D_acc = sum(this.Hessian_D_k_cell);%zeros_(1);
            count = 1;
            for index_Hessian = 1:this.Numblock
                if any(this.Hessian_u_k_cell{index_Hessian}(:))
                    V_acc(:,count) = this.Hessian_u_k_cell{index_Hessian}(:);
                    count = count+1;
                else
                    continue;
                end
            end

            if count>1
                V_acc = V_acc(:,1:count-1);
                V_inv = inv(-eye(size(V_acc,2))+(V_acc'*V_acc)/D_acc);
                % Bx_inv: vector to vector
                Bx_inv = @(x)(x(:)/D_acc-V_acc*(V_inv/(D_acc)^2)*(V_acc'*x(:)));
                sigma_max = this.Power_Iteration(Bx_inv,prod(this.sizein));
            else
                Bx_inv = @(x)(x(:)/D_acc);
                sigma_max = 1/D_acc;
            end

            if isempty(this.P)
                this.P = zeros_(this.D.sizeout);
                R = this.P;
            else
                R = this.P;
            end

            t_k_1 = 1;
            v_out = zeros_(this.D.sizein);
            if count>1
                x0 = zeros_(count-1,1);
            end
            for iterin = 1:this.maxiterin
                %%%%%%%%%%
                % Computing the gradient of the objective function
                % See Eq. (4.9) in the reference.
                v_old = v_out;
                Pold = this.P;
                v_in = reshape(Bx_inv(v-alph*(this.D'*R)),this.sizein);
                if count>1 % we have no-zeros u.
                    Eval_fun = @(x)(V_acc'*(v_in(:)-max(v_in(:)+V_acc*(x/D_acc),0))+x);
                    iter_count = 0;
                    xxx = x0;
                    %grad_norm_set = [];
                    while 1
                        iter_count = iter_count+1;
                        grad = Eval_fun(xxx);
                        temp1 =  max(v_in(:)+(V_acc*xxx)/D_acc,0);
                        index = temp1==0;
                        V_acc_temp  = V_acc; % we cannot overwrite the V_acc
                        V_acc_temp(index,:) = [];
                        Jaco_eval = -(V_acc_temp'*V_acc_temp)/D_acc+eye(count-1);
                        xxx = xxx-Jaco_eval\grad;
                        if norm(grad)<=1e-6 || iter_count<=this.MaxIterProj
                            break;
                        elseif iter_count>this.MaxIterProj
                            error('The projection reaches the maximal iteration and does not convergence.\n')
                        end
                    end
                    v_out = reshape(max(v_in(:)+V_acc*(xxx/D_acc),0),this.sizein);
                else
                    v_out = max(v_in,0);
                end
                re = norm(v_old(:)-v_out(:))/norm(v_out(:));
                if re < this.xtolTV
                    break;
                end
                %%%%%%%%%%
                % Taking a step towards minus of the gradient
                % -------------------------------------------
                % *** we may conduct line search here. ***
                % -------------------------------------------
                this.P = R+1/(this.D.norm^2*alph*sigma_max)*(this.D*(v_out));
                %%%%%%%%%%
                % Peforming the projection step
                if strcmp(this.TVtype,'iso')
                    this.P = this.P./repmat(max(1,sqrt(sum(this.P.^2, this.ndims + 1))),[ones(1,this.ndims), this.ndims]);%Project L2 ball
                elseif strcmp(this.TVtype,'l1')
                    this.P=this.P./(max(abs(this.P),1,this.ndims+1));
                else
                    error('unknown type of total variation. should be iso or l1');
                end
                %%%%%%%%%%
                % Updating R and t
                t_k = t_k_1;
                t_k_1=(1+sqrt(1+4*t_k^2))/2;
                R = this.P+(t_k-1)/(t_k_1)*(this.P-Pold);
            end
        end

        function [lambda,b] = Power_Iteration(this,Ax,size_in,tolerance)

            if nargin <=3
                tolerance = 1e-4;
            end

            b_k = randn(size_in,1);
            Ab_k = Ax(b_k);
            norm_b_k = norm(Ab_k);
            while 1
                b_k = Ab_k/norm_b_k;
                Ab_k = Ax(b_k);
                norm_b_k_1 = norm(Ab_k);
                if abs(norm_b_k_1-norm_b_k)<=tolerance
                    break;
                else
                    norm_b_k = norm_b_k_1;
                end
            end
            b = b_k;
            lambda = b'*Ab_k/(b'*b);
        end

        function v_out = ProxTVWarm(this,v,D_inv,u_inv)
            % Hessian matrix B is B = D \pm u*u'; -- test on positive first.
            % -------------------------------------------------------------------------
            alph = this.a_k*this.reg;

            if isinf(this.bounds(2)) % nonegative constraints
                if any(u_inv(:))
                    project = @(x)(reshape(this.proj_rank1_Rplus_Tao(x(:),D_inv,[],u_inv(:)),this.sizein));
                else
                    project = @(x)(reshape(this.proj_rank1_Rplus_Tao(x(:),D_inv),this.sizein));
                end
            else % box constraints
                if any(u_inv(:))
                    % weighted proximal operator HERE
                    % Prox_h^B(w(p,q))
                    project = @(x)(reshape(this.proj_rank1_box_Tao(x(:),D_inv,[],u_inv(:)),this.sizein));
                else
                    project = @(x)(reshape(this.proj_rank1_box_Tao(x(:),D_inv),this.sizein));
                end
            end
            if any(u_inv(:))
                sigma_max = D_inv+sum(u_inv.*u_inv,'all');
            else
                sigma_max = D_inv;
            end
            % One can set larger stepsize - Ref. Results of my new paper. or set 3/4 to be 1 as we used before.
            if any(u_inv(:))
                % formulate operator Bx or B^{-1}x
                Binvx = @(x)(D_inv*x+sum(u_inv.*x,'all')*u_inv);
            else
                Binvx = @(x)(D_inv*x);
            end

            if isempty(this.P)
                this.P = zeros_(this.D.sizeout);
                R = this.P;
            else
                R = this.P;
            end

            t_k_1 = 1;
            v_out = zeros_(this.D.sizein);
            %fun_set = [];
            for iterin = 1:this.maxiterin
                v_old = v_out;
                Pold = this.P;
                % weighted proximal operator HERE
                % Prox_h^B(w(p,q))
                v_in = v-alph*Binvx(this.D'*R);
                v_out = project(v_in);
                %                 end
                re = norm(v_old(:)-v_out(:))/norm(v_out(:));
                if re < this.xtolTV
                    break;
                end
                %%%%%%%%%%
                % Taking a step towards minus of the gradient
                % -------------------------------------------
                % *** we may conduct line search here. ***
                % -------------------------------------------
                this.P = R+1/(this.D.norm^2*alph*sigma_max)*(this.D*(v_out));

                %%%%%%%%%%
                % Peforming the projection step
                if strcmp(this.TVtype,'iso')
                    this.P = this.P./repmat(max(1,sqrt(sum(this.P.^2, this.ndims + 1))),[ones(1,this.ndims), this.ndims]);%Project L2 ball
                elseif strcmp(this.TVtype,'l1')
                    this.P=this.P./(max(abs(this.P),1,this.ndims+1));
                else
                    error('unknown type of total variation. should be iso or l1');
                end
                %%%%%%%%%%
                % Updating R and t
                t_k = t_k_1;
                t_k_1=(1+sqrt(1+4*t_k^2))/2;

                R = this.P+(t_k-1)/(t_k_1)*(this.P-Pold);

            end
        end

        function varargout = proj_rank1_box_Tao(this,varargin)
            % PROJ_RANK1_BOX returns the scaled proximity operator for box constraints
            prox            = @(x,t) max( min(this.bounds(2),x), this.bounds(1));
            prox_brk_pts    = @(s) [this.bounds(1),this.bounds(2)]; % since projection, scaling has no effect
            [varargout{1:nargout}] = this.prox_rank1_generic_Tao(prox, prox_brk_pts, varargin{:} );
        end

        function varargout = proj_rank1_Rplus_Tao(this,varargin)
            prox            = @(x,t) max(0, x);
            prox_brk_pts    = @(s) 0;
            [varargout{1:nargout}] = this.prox_rank1_generic_Tao(prox,prox_brk_pts,varargin{:} );
        end

        function [x,a,cnt] = prox_rank1_generic_Tao(this,prox,prox_brk_pts,x0,D,u,uinv,lambda,linTerm, plusminus,INVERT)
            n   = length(x0);

            if nargin < 6 || isempty(u), u = 0; end
            if nargin<7 || isempty(u), uinv = 0; end
            if nargin < 8, lambda = []; end
            if nargin < 9, linTerm = []; end
            if nargin < 10 || isempty(plusminus), plusminus = 1; end
            assert( plusminus==-1 | plusminus==+1 )
            if nargin < 11 || isempty(INVERT), INVERT = true; end

            if size(D,2) > 1, d = diag(D); else d = D; end % extract diagonal part
            if any( d < 0 ), error('D must only have strictly positive entries'); end

            if all( u==0 )
                % Just a diagonal scaling, so this code is overkill,
                % but we should be able to handle it for consistency.
                NO_U = true;
            else
                NO_U = false;
                if numel(u) > length(u)
                    error('u must be a vector, not a matrix');
                end
            end

            if plusminus < 0 && all( d==d(1) )
                minE = d(1) + plusminus*norm(u)^2;
                if minE <= 0, error('The scaling matrix is not positive definite'); end
            end

            % In all cases, we find prox_h^V, but how we define V
            %   in terms of d and u depends on "INVERT"
            if INVERT
                % So V^{-1} = diag(d)     + sigma*u*u'
                % and     V = diag(1./d)  - sigma*uinv*uinv';
                Vinv = @(y) d.*y + (plusminus*(u'*y))*u;

                %   The code below expects V = diag(dd) + sigma*uu*uu', so...
                dd          = 1./d;
                uu          = uinv;
                plusminus   = -plusminus;

                % The code also requires uu./dd and 1./dd, so define these here
                ud      = uu./dd;
                %ud      = u/sqrt(1+u'*(u./d)); % more accurate? % 6.01e-3 error
                dInv    = 1./dd;
            else
                % Here, V    = diag(d) + sigma*u*u'
                % and V^{-1} = diag(1./d) - sigma*uinv*uinv';
                Vinv = @(y) y./d - (plusminus*(uinv'*y))*uinv;

                %   The code below expects V = diag(dd) + sigma*uu*uu', so...
                dd          = d;
                uu          = u;
                %plusminus   = plusminus;

                % The code also requires uu./dd and 1./dd, so define these here
                ud      = uu./dd;
                dInv    = 1./dd;
            end
            if NO_U, uu = 0; ud = 0; end % any value, since we won't use them...
            if ~isempty(lambda)

                if any(lambda==0), error('scaling factor lambda must be non-zero'); end
                % note that lambda < 0 should be OK
                x0 = lambda.*x0;

                % Scale V = diag(dd) + sigma*uu*uu' by V <-- diag(1./lambda)*V*diag(1./lambda)
                dd = dd./(lambda.^2);
                uu = uu./lambda;
                ud = ud.*lambda;
                dInv    = 1./dd;
            end
            t   = prox_brk_pts(1./dd);
            if size(t,1) < n
                if size(t,1) > 1
                    error('"prox_brk_pts" should return a ROW VECTOR of break points');
                end
                % otherwise, assume each component identical, so scale
                t = repmat(t,n,1);
            end
            if ~isempty(linTerm) && norm(linTerm)>=0
                if isempty(lambda)
                    x0  = x0 - Vinv(linTerm);
                else
                    % V is scaled V <-- diag(1./lambda)*V*diag(1./lambda)
                    %   so Vinv is scaled the opposite.
                    % linTerm is scaled linTerm <== linTerm./lambda
                    x0  = x0 - lambda.*Vinv(linTerm);
                end
            end

            % The main heart:
            X       = @(a) prox( x0 - (plusminus*a)*ud, dInv );

            % Early return if we have only a diagonal scaling...
            if NO_U
                % in this case, "alpha" is irrelevant
                x   = prox( x0, dInv );
                if ~isempty(lambda)
                    % Undo the scaling of x <-- lambda.*x
                    x = x./lambda;
                end
                return;
            end

            brk_pts = bsxfun(@times,((plusminus*dd)./uu),  bsxfun(@minus,x0,t) );
            % we use nonegative that we do not need this.
            brk_pts = brk_pts(~isinf(brk_pts)); % in case lwr/upr=Inf for box
            brk_pts = unique(brk_pts(:)); % will sort and remove duplicates
            

            % Main for-loop:
            % "lower bound" are "a" for which p <= 0
            % "upper bound" are "a" for which p >= 0
            % if a is increasing, so is p(a) (double-check for both plusminus cases )
            lwrBnd       = 0;
            uprBnd       = length(brk_pts) + 1;
            iMax         = ceil( log2(length(brk_pts)) ) + 1;

            for i = 1:iMax
                if uprBnd - lwrBnd <= 1
                    %dispp('Bounds are too close; breaking');
                    break;
                end
                j = round(0.5*(lwrBnd+uprBnd));
                %printf('j is %d (bounds were [%d,%d])\n', j, lwrBnd,uprBnd );
                if j==lwrBnd
                    %dispp('j==lwrBnd, so increasing');
                    j = j+1;
                elseif j==uprBnd
                    %dispp('j==uprBnd, so increasing');
                    j = j-1;
                end

                % to decide which a let p = 0;
                a   = brk_pts(j);
                x   = X(a);  % the prox
                p   = a + dot(uu,x0-x);

                if p > 0 % search one a let p = 0;
                    uprBnd = j;
                elseif p < 0
                    lwrBnd = j;
                end
            end
            cnt     = i; % number of iterations we took
            % Now, determine linear part, which we infer from two points.
            % If lwr/upr bounds are infinite, we take special care
            % e.g., we make a new "a" slightly lower/bigger, and use this
            % to extract linear part.
            if lwrBnd == 0
                a2 = brk_pts( uprBnd );
                a1 = a2 - 10; % arbitrary
                aBounds = [-Inf,a2];
            elseif uprBnd == length(brk_pts) + 1
                a1 = brk_pts( lwrBnd );
                a2 = a1 + 10; % arbitrary
                aBounds = [a1,Inf];
            else
                % In general case, we can infer linear part from the two break points
                a1 = brk_pts( lwrBnd );
                a2 = brk_pts( uprBnd );
                aBounds = [a1,a2];
            end
            x1 = X(a1);
            x2 = X(a2);
            dx = (x2 - x1)/(a2-a1);
            % Thus for a in (a1,a2), x(a) = x1 + (a-a1)*dx
            % Solve 0 = a + dot( uu, y - (x1 + (a-a1)*dx ) )
            %         = a + dot(uu,y - x1 + a1*dx ) - a*dot(uu,dx)
            % so:
            a = dot( uu, x0 - x1 + a1*dx)/( -1 + dot(uu,dx) );
            if a < aBounds(1) || a > aBounds(2), error('alpha is not in the correct range'); end
            % If we were not linear, we could do a root-finding algorithm, e.g.,
            % a = fzero( @(a) a+dot(uu,x0-X(a)), a );

            % Now, the solution:
            x = X(a);

            if ~isempty(lambda)
                x = x./lambda;
            end
        end
    end
end
