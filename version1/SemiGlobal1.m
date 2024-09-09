function [U, mniter, matvecs, est_errors, history] = SemiGlobal1(Gop, Gdiff_op, Gdiff_matvecs,...
    ihfun, ev_domain, ui, tgrid, Nts, Nt_ts, Nfm, tol, options, data, varargin)
% The program solves time-dependent Schroedinger equation for a
% time-dependent, nonlinear Hamiltonian, using the semi-global propagator 
% by Hillel Tal-Ezer.
% There is an option to include an inhomogeneos term.
%
% Input:
% Gop: A function handle of the form @(u, t, v, more arguments). Returns 
% G(u, t)v (in quantum mechanics, G(u, t) = -iH(u, t)/hbar). The "more arguments"
% will be written in the place of "varargin" in this program.
% Gdiff_op: A function handle of the form @(u1, t1, u2, t2, more arguments). Returns
% (G(u1, t1) - G(u2, t2))u1. u1 and t1 represent several time-points in 
% separated columns, while u2 and t2 represent a single time-point.
% Note: "more arguments" must be the same as for Gop.
% Gdiff_matvecs: The number of matrix-vector multiplications required for
% the computation of the operation of Gdiff_op for each time-point (usually less than 2).
% (The number of matrix-vector multiplications is counted as the number of
% large scale operations with the highest scaling with the dimension of the
% problem.)
% ihfun: A function handle of the form: @(t, more arguments). Returns the
% inhomogeneous term vector. If there is no inhomogeneous term, put []
% instead.
% Note: "more arguments" must be the same as for Gop.
% ev_domain: The boundaries of the eigenvalue domain of G, when a Chebyshev
% algorithm is used for the computation of the function of matrix. For an
% Arnoldi algorithm, put [].
% ui: the initial state vector.
% tgrid: the time grid of the desired output solution, U. Has to be ordered in an
% increasing order for a forward propagation, and a decreasing order for a
% backward propagation. Has to contain at least two time-points: The
% initial and final point.
% Nts: the number of time steps of the propagation. 
% Nt_ts: the number of interior Chebyshev time points in each time-step, used during
% the computational process (M in the article).
% Nfm: the number of expansion terms for the computation of the function of
% matrix (K in the article).
% tol: the desired tolerance of the total propagation error; for a manual
% mode with a predefined iteration number, insert [] instead. In this case,
% The iteration number of the first time-step is given by options.Niter1st,
% and the iteration number for all other time-steps is given by
% options.Niter.
% options: A structure of propagation options. The function SGdefault_op.m
% sets the defaults. See the structure fields therein.
% data: The required data computed as a preparation stage. Can be produced by the
% function SGdata.m.
% If the values of options and/or data are not inserted, they are computed
% as the default in the program.
% Note: If the propagation itself consumes a small amount of time, but the
% overall computation consists of many propagations (as in optimal control
% problems), it's not recommended not to insert the variables options and data as an input.
%
% Output:
% U: contains the solution at the time points specified in tgrid, in separate columns.
% mniter: the mean number of iteration for a time step, excluding the first step,
% which typically requires more iterations. Usually, should be 1 for ideal
% efficiency of the algorithm.
% matvecs: Structure of the number of Hamiltonian operations. It has the
% following fields:
%   propagation: The number of Hamiltonian operations required for the
%   propagation process;
%   error: The number of Hamiltonian operations required for the error
%   estimation;
%   total: The total number of Hamiltonian operations (the sum of the two
%   other fields).
% est_errors: A structure that contains the total estimated relative errors, based 
% on the assumption of additive error accumulation during the propagtion.
% It has the following fields:
%   texp_cheap: The estimated relative error resulting from the time-expansions in
%   each time-step; the estimation is relatively cheap numerically, but
%   overestimates the error, typically by 2-3 orders of magnitude.
%   texp_exact: The estimated relative error resulting from the time-expansions in
%   each time-step; the resulting estimation is much more precise than 
%   texp_cheap. However, in the case of odd Nt_ts the estimation is more
%   expansive numerically. Not computed for odd Nt_ts if options.texp_er_odd = false
%   (the default is options.texp_er_odd = true).
%   fm: The estimated relative error resulting from the
%   computation of the function of matrix.
%   conv: The total estimated relative convergence error. It is computed
%   based on the different convergence error estimations
%   for each time-step. The sums of the different estimations are represented
%   by the following fields.
%   conv_cheb: The total estimated relative convergence error based on
%   numerical integration at the Chebyshev points. Becomes demanding for
%   Gdiff_matvecs>0. The default is not to perform this computation in this
%   case (see SGdefault_op.m).
%   conv_texp: Based on multiplication of the integration by the end-point
%   of the interval by a correction factor; the correction factor is based
%   on the assumption that the extrapolation error is dominated by the
%   time-expansion error.
%   conv_fm: Based on multiplication of the integration by the end-point
%   of the interval by a correction factor; the correction factor is based
%   on the assumption that the extrapolation error is dominated by the
%   function of matrix computation error.
%   total: The total estimated error; it is the sum of the fields
%   texp_exact, fm and conv. In case that texp_exact isn't computed, it is
%   replaced by texp_cheap.
% history: A structure which contains the history of propagation with the
% following fields:
%   t: The time grid of propagation
%   U: The solution at history.t
%   texp_error_cheap: The cheap time-expansion error estimation for all time-steps
%   texp_error_exact: The exact time-expansion error estimation for all time-steps
%   fm_error: The function of matrix error estimation for all time-steps
%   conv_error, conv_error_cheb, conv_error_texp, conv_error_fm: The convergence error
%   estimations for all time-steps (see the description for the fields of
%   the output structure est_errors).
%   total_error: The total error estimation for all time-steps (as
%   described in the corresponding field of est_errors)
%   reldif: The relative difference of the solution from the previous 
%   iterated solution for all time-steps. Provides an estimation of the
%   error of the previous iterated solution.
%   niter: The number of iterations for each time-step
%   stab_factor: The stability factor; it is a scalar for the Chebyshev
%   algorithm. For the Arnoldi algorithm it is computed for each time-step.
%   Problematic values are in the order of ~0.1-1 or higher.
% Called functions: SGdefault_op, SGdata, chebcM, chebc2result, vchebMop,
% chebweights, createKrop, getRvKr, divdif, new_divdif

% Author: Ido Schaefer
% ido.schaefer@gmail.com

%   Setting defaults:
    if nargin<12
        options = SGdefault_op;
    end
    if nargin<13
        data = SGdata(options);
    end
    if ~isempty(tol)
        % The tolerance is specified by the user; allowed error per time-step:
        tol_ts = tol/Nts;
        % The number of iterations is determined adaptively by the allowed tolerance:
        tol_mode = true;
    else
        % Employing a manual mode for determination of the number of
        % iterations. The number of iterations is predefined,
        % options.Niter1st for the 1st time-step, and options.Niter for
        % the rest of the time-steps.
        tol_mode = false;
    end
    if isempty(options.tol_f)
        % The default tolerance of the f_{Nt_ts}(z, t) computation:
        if tol_mode
            options.tol_f = tol_ts;
        else
            options.tol_f = eps;
        end
    end
    if isempty(options.tol1st)
        tol1st_mode = false;
    else
        % The 1st time-step is treated differently than the rest of the
        % time-steps with its own tolerance parameter:
        tol1st_mode = true;
    end
    if isempty(options.conv_er_cheb)
        % The default for the computation of the convergence error by
        % integration at Chebyshev points; the computation won't be
        % performed if the evaluation of the extended "inhomogeneous" vector costs
        % non-negligible computational effort:
        if Gdiff_matvecs == 0
            options.conv_er_cheb = true;
        else
            options.conv_er_cheb = false;
        end
    end
    % If the eigenvalue domain is not specified, the Arnoldi approach is
    % employed.
    Arnoldi = isempty(ev_domain);
    % In order to detect if the propagation is a forward or a backward
    % propagation:
    direction = sign(tgrid(2) - tgrid(1));
    Nt = length(tgrid);
    tinit = tgrid(1);
    tf = tgrid(Nt);
    % The length of the time interval of the whole propagation (can be negative):
    T = tf - tinit;
    % The length of the time step interval:
    Tts = T/Nts;
    Nu = length(ui);
    U = zeros(Nu, Nt);
    U(:, 1) = ui;
    % The Chebyshev points for expansion in time, in the domain in which the Chebyshev expansion
    % is defined: [-1 1]
    tcheb = -cos(((1:Nt_ts) - 1)*pi/(Nt_ts-1));
    % The Chebyshev points for expansion in time in the domain of the time
    % variable:
    t_ts = 0.5*(tcheb + 1)*Tts;
    if options.conv_er_cheb
        % Chebyshev integration weights for caculation of est_errors.conv_cheb:
        wcheb = chebweights(Nt_ts, Tts).';
    end
    % The parity of Nt_ts:
    Nt_ts_is_even = (mod(Nt_ts, 2) == 0);
    if Nt_ts_is_even
        % The index of the time-point in the middle of the time step; for
        % even Nt_ts it is the test-point:
        tmidi = Nt_ts + 1;
        % For time-expansion error computation:
        Etexp_factor_even = 8*abs(Tts/((Nt_ts^2 - 1)*(Nt_ts - 3)));
        % s_ext evaluation time points for error estimation (just one point
        % for even Nt_ts):
        texp_er_tpoints = Tts/2;
        % The index of the test-point:
        texp_er_i = Nt_ts + 1;
        % The indices for the computation of \bar{G}(t)u(t) in s_ext (without the test-points;
        % note that this differs from SemiGlobal.m):
        s_ext_i = 1:Nt_ts;
    else
        % For odd Nt_ts, the middle time-point is also in the middle of the
        % time-step:
        tmidi = ceil(Nt_ts/2);
        if options.texp_er_odd
            % For the exact time-expansion error computation; performed
            % only if specified by the field texp_er_odd in options:
            Etexp_factor_odd = 4*abs(Tts*(Nt_ts - 1)/(Nt_ts*(Nt_ts^2 - 4)*(Nt_ts - 4)));
            % The error estimaiton for odd Nt_ts requires two test points:
            texp_er_tpoints = [(t_ts(tmidi) + t_ts(tmidi + 1))/2, (t_ts(tmidi) + t_ts(tmidi - 1))/2];
            % The indices of the test-points:
            texp_er_i = [Nt_ts + 1, Nt_ts + 2];
        else
            % Only the cheap time-expansion error estimation is performed: 
            texp_er_tpoints = (t_ts(tmidi) + t_ts(tmidi + 1))/2;
            texp_er_i = Nt_ts + 1;
        end
        % If Nt_ts is odd, the midpoint needn't be computed:
        s_ext_i = [1:(tmidi - 1), (tmidi + 1):Nt_ts];
        if options.conv_er_cheb
            % For odd Nt_ts, the value of the integrand at the middle point is
            % zero, and thus not included in the Chebyshev weights.
            wcheb(tmidi) = [];
        end
    end
    % For the cheap time-expansion error estimation:
    Etexp_factor_cheap = 4*Tts;
    % The total number of internal time-points, including test-points:
    total_Nt_ts = Nt_ts + length(texp_er_tpoints);
    % The number of points in s_ext_i:
    Ns_ext_i = length(s_ext_i);
    % The interior time points of the current time step for interpolation,
    % and the next time step for extrapolation of the guess solution into the next step:
    t_2ts = [t_ts, texp_er_tpoints, Tts + t_ts(2:Nt_ts), Tts + texp_er_tpoints];    
    % The full propagation grid:
    propagation_grid = [kron((tinit:Tts:(tf - Tts)), ones(1, total_Nt_ts - 1)) + ...
        kron(ones(1, Nts), [t_ts(1:(Nt_ts - 1)), texp_er_tpoints]), tf];
    if nargout>4
        if ~options.save_memory
            history.U = zeros(Nu, Nts*(Nt_ts - 1) + 1);
            history.U(:, 1) = ui;
        end
        history.t = [kron(tinit:Tts:(T - Tts), ones(1, Nt_ts - 1)), tf] + [kron(ones(1, Nts), t_ts(1:(Nt_ts - 1))), 0];
        history.texp_error_cheap = zeros(1, Nts);
        if Nt_ts_is_even || options.texp_er_odd
            history.texp_error_exact = zeros(1, Nts);
        end
        %history.f_error = zeros(1, Nts);
        history.fm_error = zeros(1, Nts);
        history.reldif = zeros(1, Nts);
        history.conv_error = zeros(1, Nts);
        history.conv_error_texp = zeros(1, Nts);
        history.conv_error_fm = zeros(1, Nts);
        if options.conv_er_cheb
            history.conv_error_cheb = zeros(1, Nts);
        end
        history.niter = zeros(1, Nts);
    end
    % Necessary for error estimation of f_{Nt_ts}(z, t):
    factorialNt_ts = factorial(Nt_ts);
    if ~Arnoldi
        % If the eigenvalue domain is specified, a Chebyshev approximation
        % for the function of matrix is employed.
        min_ev = ev_domain(1);
        max_ev = ev_domain(2);
        % Computing the coefficients for the Chebyshev expansion of the
        % function of matrix, in all the interior time points.
        % CchebFts contains the coefficients of the current time step, and
        % CchebFnext contains the coefficients of the next one.
        Ccheb_f_comp = f_chebC(t_2ts(2:(2*total_Nt_ts - 1)));
        Ccheb_f_ts = Ccheb_f_comp(:, 1:(total_Nt_ts - 1));
        Ccheb_f_next = Ccheb_f_comp(:, total_Nt_ts:(2*total_Nt_ts - 2));   
        % Computing sampling points for the error test of the
        % function of matrix error at the maxima of the Chebyshev polynomial:
        cheb_testp = -cos((0:Nfm)*pi/Nfm).';
        ztest = (cheb_testp*(max_ev - min_ev) + min_ev + max_ev)/2;
        fztest = f_fun(ztest, Tts);
        f_scalar_error = max(abs(chebc2result(Ccheb_f_ts(:, Nt_ts - 1), ev_domain, ztest) - fztest));
        % Stability factor:
        stability_factor = f_scalar_error*max(abs(ev_domain))^Nt_ts/factorialNt_ts;
        % Estimated stability criterion (a boolean variable):
        instability = stability_factor>0.1;
        if instability
            fprintf('Warning: Instability in the propagation process may occur.\n')
        end
        if nargout>4
            history.stab_factor = stability_factor;
        end
    else
        % For the Arnoldi algorithm, the stability test is performed in
        % each time step.
        instability = false;
        if nargout>4
            history.stab_factor = zeros(1, Nts);
        end
    end
    % A boolean that indicates if the convergence has failed in at least one
    % time-step:
    conv_failure = false;
    est_errors.texp_cheap = 0;
    if Nt_ts_is_even || options.texp_er_odd % Condition for computation of the exact estimation
        est_errors.texp_exact = 0;
    end
    est_errors.fm = 0;
    est_errors.conv = 0;
    est_errors.conv_texp = 0;
    est_errors.conv_fm = 0;
    if options.conv_er_cheb
        est_errors.conv_cheb = 0;
    end
    % Computing the matrix of the time Taylor polynomials.
    % timeMts contains the points in the current time step, and timeMnext
    % contains the points in the next time step:
    timeMcomp = maketM(t_2ts(2:(2*total_Nt_ts - 1)), Nt_ts);
    timeMts = timeMcomp(:, 1:(total_Nt_ts - 1));
    timeMnext = timeMcomp(:, total_Nt_ts:(2*total_Nt_ts - 2));    
    % Computing the coefficients of the transformation from the Newton
    % interpolation polynomial terms, to a Taylor like form:
    Cr2t = r2Taylor4(t_ts, Tts);
    % The extended "inhomogeneous" vectors:
    s_ext = zeros(Nu, total_Nt_ts);
    % Newton coefficients for the time-expansion of s_ext:
    Cnewton = zeros(Nu, total_Nt_ts);
    % The v vectors are defined recursively, and contain information about
    % the time dependence of the s_ext vectors:
    v_vecs = zeros(Nu, Nt_ts + 1);
    % If there is no inhomogeneous term in the equation, ihfun is empty, and there_is_ih = false.
    % If there is, there_is_ih == true.
    there_is_ih = ~isempty(ihfun);
    if there_is_ih
        s = zeros(Nu, total_Nt_ts);
        s(:, 1) = ihfun(tinit, varargin{:});
    end
    % The 0'th order approximation is the first guess for the first time step.
    % Each column represents an interior time point in the time step:
    Uguess = guess0(ui, total_Nt_ts);
    Unew = zeros(Nu, total_Nt_ts);
    % The total number of iterations, excluding the first time step, for
    % calculation of mniter:
    allniter = 0;
    % These variables are used to determine which points in tgrid are in the computed time-step.  
    tgrid_lowi = 2;
    tgrid_upi = 1;
    for tsi = 1:Nts
        % The time of the interior time points within the time-step:
        t = propagation_grid((tsi - 1)*(total_Nt_ts - 1) + [(1:(Nt_ts - 1)), total_Nt_ts, texp_er_i - 1]);
        % The first guess for the iterative process, for the convergence of the u
        % values. Each column represents an interior time point in the time step:
        Ulast = Uguess;
        Unew(:, 1) = Ulast(:, 1);
        v_vecs(:, 1) = Ulast(:, 1);
        if there_is_ih
            % Computing the inhomogeneous term:
            s(:, 2:total_Nt_ts) = ihfun(t(2:total_Nt_ts), varargin{:});
        end
        if options.conv_er_cheb
            % Computation of \bar{G}u at the Chebyshev points for
            % integration. Note that for odd Nt_ts, \bar{G}u = 0, and therefore not
            % calculated.
            Gdiff_u_new = Gdiff_op(Ulast(:, s_ext_i), t(s_ext_i), Ulast(:, tmidi), t(tmidi), varargin{:});
        else
            % Computation of \bar{G}u at the end-point of the time-step
            % only:
            Gdiff_u_new_end = Gdiff_op(Ulast(:, Nt_ts), t(Nt_ts), Ulast(:, tmidi), t(tmidi), varargin{:});
        end
        % Starting an iterative process until convergence:
        niter = 0;
        conv_error = Inf;
        while (tsi>1 && niter<options.Niter && ((tol_mode && conv_error>tol_ts) || ~tol_mode)) ||...
              (tsi == 1 && niter<options.Niter1st &&...
                ((~tol1st_mode && ((tol_mode && conv_error>tol_ts) || ~tol_mode)) || (tol1st_mode && conv_error>options.tol1st)))
            % Setting the inhomogeneous s_ext vectors. 
            if options.conv_er_cheb
                s_ext(:, s_ext_i) = Gdiff_u_new;
            else
                s_ext(:, s_ext_i) =...
                    [Gdiff_op(Ulast(:, s_ext_i(1:(Ns_ext_i - 1))), t(s_ext_i(1:(Ns_ext_i - 1))), Ulast(:, tmidi), t(tmidi), varargin{:}),...
                    Gdiff_u_new_end];
                % Note that for odd Nt_ts, s_ext(:, tmidi) is equivalent to
                % s(:, tmidi), and therefore not calculated.
            end    
            if there_is_ih
                % Saving the previous \bar{G}u for calculation of the
                % convergence error. Required only when s_ext is different
                % than \bar{G}u - when there's a source term.
                if options.conv_er_cheb
                    Gdiff_u = Gdiff_u_new;
                else
                    Gdiff_u_end = Gdiff_u_new_end;
                end
                if ~Nt_ts_is_even
                    % For odd Nt_ts, the extended "source term" at the middle
                    % time-point is just the source term:
                    s_ext(:, tmidi) = s(:, tmidi);
                end
                % The inhomogeneous term is added to the s_ext vectors:
                s_ext(:, s_ext_i) = s_ext(:, s_ext_i) + s(:, s_ext_i);
            end
            % Calculation of the coefficients of the form of Taylor
            % expansion, from the coefficients of the Newton
            % interpolation at the points t_ts.
            % The divided differences are computed by the function divdif.
            % For numerical stability, we have to transform the time points
            % in the time step, to points in an interval of length 4:
            [Cnewton(:, 1:Nt_ts), diagonal] = divdif(t_2ts(1:Nt_ts)*4/Tts, s_ext(:, 1:Nt_ts));
            % Calculating the Taylor like coefficients:
            Ctaylor = zeros(Nu, Nt_ts);
            Ctaylor(:, 1) = Cnewton(:, 1);
            for Newtoni = 2:Nt_ts
                Ctaylor(:, 1:Newtoni) = Ctaylor(:, 1:Newtoni) + Cnewton(:, Newtoni)*Cr2t(Newtoni, 1:Newtoni);
            end
            % Calculation of the v vectors:
            for polyi = 2:(Nt_ts+1)
                v_vecs(:, polyi) = (Gop(Ulast(:, tmidi), t(tmidi), v_vecs(:, polyi-1), varargin{:}) + Ctaylor(:, polyi-1))/(polyi - 1);
            end
            if ~min(isfinite(v_vecs(:, Nt_ts + 1)))
                % It means that the algorithm diverges.
                % In such a case, change Nts, Nt_ts and/or Nfm.
                fprintf('\nError: The algorithm diverges (in time step No. %d).\n', tsi);
                mniter = allniter/(tsi - 1);
                if tsi == 1
                    matvecs = [];
                end
                return
            end
            if Arnoldi
                % Creating the Krylov space by the Arnoldi iteration procedure,
                % in order to approximate f(G,t)v_vecs(:, Nt_ts + 1):
                [Upsilon, Hessenberg] = createKrop(@(v) Gop(Ulast(:, tmidi), t(tmidi), v, varargin{:}), v_vecs(:, Nt_ts + 1), Nfm);
                % Obtaining eigenvalues of the Hessenberg matrix:
                eigval = eig(Hessenberg(1:Nfm, 1:Nfm));
                % The test point is the average point of the eigenvalues:
                avgp = sum(eigval)/Nfm;
                samplingp = [eigval; avgp];
                capacity = get_capacity(eigval, avgp);
                % Obtaining the expansion vectors for a Newton approximation of f_{Nt_ts}(G,t)v_vecs(:, Nt_ts + 1)
                % in the reduced Krylov space:
                RvKr = getRv(v_vecs, Hessenberg, samplingp, capacity);
                % Calculation of the solution at all time points
                % within the time step:
                [Unew(:, 2:total_Nt_ts), f_abs_error, fUerror] = Ufrom_vArnoldi(v_vecs, timeMts, Upsilon, RvKr, samplingp, capacity);
                %[Unew(:, 2:total_Nt_ts), f_error, fUerror] = Ufrom_vArnoldi(v_vecs, timeMts, Upsilon, RvKr, samplingp, capacity);
            else
                % Employing a Chebyshev approximation for the function of
                % matrix computation.
                % Vcheb is a matrix. Contains the vectors:
                % T_n(G(Ulast(:,tmidi), t(tmidi)))*v_vecs(: ,Nt_ts + 1),  n = 0, 1, ..., Nfm-1
                % where the T_n(z)'s are the Chebyshev polynomials.
                % The n'th vector is the (n+1)'th column of Vcheb.
                Vcheb = vchebMop(@(v) Gop(Ulast(:, tmidi), t(tmidi), v, varargin{:}), v_vecs(:, Nt_ts + 1), min_ev, max_ev, Nfm);
                % Calculation of the solution at all time points
                % within the time step:
                [Unew(:, 2:total_Nt_ts), fUerror] = Ufrom_vCheb(v_vecs, timeMts, Vcheb, Ccheb_f_ts);
                %[Unew(:, 2:total_Nt_ts), f_error, fUerror] = Ufrom_vCheb(v_vecs, timeMts, Vcheb, Ccheb_f_ts);
            end
            % The relative difference from the previous solution in the
            % iterative process:
            reldif = norm(Unew(:, Nt_ts) - Ulast(:, Nt_ts))/norm(Ulast(:, Nt_ts));
            if niter == 0
                reldif1st = reldif;
            end
            % Convergence error estimation:
            if options.conv_er_cheb
                % Calculation of the new \bar{G}u vectors:
                Gdiff_u_new = Gdiff_op(Unew(:, s_ext_i), t(s_ext_i), Unew(:, tmidi), t(tmidi), varargin{:});
                if there_is_ih
                    % The cheap estimation is based on the end-point of the
                    % integration interval:
                    conv_error_cheap = norm(Gdiff_u_new(:, Ns_ext_i) - Gdiff_u(:, Ns_ext_i))*Tts/norm(Unew(:, Nt_ts));
                    % An accurate estimation based on integration at the
                    % Chebyshev points:
                    conv_error_cheb = norm((Gdiff_u_new - Gdiff_u)*wcheb)/norm(Unew(:, Nt_ts));
                else
                    % In the homogeneous case, s_ext is equivalent to the
                    % previous \bar{G}u, which needn't be stored in memory:
                    conv_error_cheap = norm(Gdiff_u_new(:, Ns_ext_i) - s_ext(:, Nt_ts))*Tts/norm(Unew(:, Nt_ts));
                    conv_error_cheb = norm((Gdiff_u_new - s_ext(:, s_ext_i))*wcheb)/norm(Unew(:, Nt_ts));
                end
                if tsi > 1
                    % Computation of the convergence error based on
                    % multiplication of the cheap estimation by correction
                    % coefficients:
                    conv_error_texp = conv_error_cheap*data.correct_texp(niter + 1, Nt_ts - 1);
                    conv_error_fm = conv_error_cheap*data.correct_fm(niter + 1, Nt_ts + Nfm - 1);
                    % The accepted convergence error is the maximum of the
                    % error by Chebyshev integration and the error based on
                    % multiplication by the correction coefficient. It is
                    % assumed that the correction coefficient that yields
                    % the closest estimation to conv_error_cheb is the
                    % right one.
                    if abs(conv_error_texp - conv_error_cheb)<=abs(conv_error_fm - conv_error_cheb)
                        % It is assumed that the extrapolation error is
                        % dominated by the time-expansion error:
                        conv_error = max([conv_error_texp, conv_error_cheb]);
                    else
                        % It is assumed that the extrapolation error is
                        % dominated by the function of matrix error:
                        conv_error = max([conv_error_fm, conv_error_cheb]);
                    end
                    % Note: Actually, it's possible to estimate both the
                    % time-expansion extrapolation error and the function
                    % of matrix extrapolation error. Then, it's possible to
                    % determine the right correction coefficient. I have no
                    % time right now for that. Hopefully in a future
                    % version.
                else
                    % Estimation for the first time-step based on
                    % correction of the cheap estimation. It isn't related
                    % to extrapolation error, and thus conv_error_texp
                    % and conv_error_fm are set to the same value:
                    conv_error_texp = conv_error_cheap*data.correct_1st(niter + 1);
                    conv_error_fm = conv_error_texp;
                    conv_error = max([conv_error_texp, conv_error_cheb]);
                end
            else
                % The convergence error by Chebyshev integration isn't
                % computed. The estimation is solely based on the cheap
                % estimation multiplied by a correction coefficient.
                Gdiff_u_new_end = Gdiff_op(Unew(:, Nt_ts), t(Nt_ts), Unew(:, tmidi), t(tmidi), varargin{:});
                if there_is_ih
                    conv_error_cheap = norm(Gdiff_u_new_end - Gdiff_u_end)*Tts/norm(Unew(:, Nt_ts));
                else
                    conv_error_cheap = norm(Gdiff_u_new_end - s_ext(:, Nt_ts))*Tts/norm(Unew(:, Nt_ts));
                end
                if tsi > 1
                    conv_error_texp = conv_error_cheap*data.correct_texp(niter + 1, Nt_ts - 1);
                    conv_error_fm = conv_error_cheap*data.correct_fm(niter + 1, Nt_ts + Nfm - 1);
                    conv_error = max([conv_error_texp, conv_error_fm]);
                else
                    conv_error_texp = conv_error_cheap*data.correct_1st(niter + 1);
                    conv_error_fm = conv_error_texp;
                    conv_error = conv_error_texp;
                end
            end
            % The solution before the last one is stored for
            % computation of the time-expansion error after the end of the
            % iterative process:
            Ulast2 = Ulast;
            Ulast = Unew;
            niter = niter + 1;
        end
        % Test for stability:
        %fprintf('%d\n', norm(v_vecs(:, Nt_ts + 1)))
        if Arnoldi && (~instability || nargout>4)
            f_scalar_error = f_abs_error/norm(v_vecs(:, Nt_ts + 1));
            stability_factor = f_scalar_error*max(abs(eigval))^Nt_ts/factorialNt_ts;
            % If a possibility of instability hasn't been detected yet, the
            % stability criterion is checked:
            if ~instability
                % Estimated stability criterion (a boolean variable):
                instability = stability_factor>0.1;
                if instability
                    fprintf('Warning: Instability in the propagation process may occur (detected in time step No. %d).\n', tsi)
                end
            end
            if nargout>4
                history.stab_factor(tsi) = stability_factor;
            end
        end
        % Detection of the first appearance of convergence failure:
        if niter>1 && reldif>reldif1st && ~conv_failure
            conv_failure = true;
            fprintf('Warning: Convergence failure (first occured in time step No. %d).\n', tsi)
            % In such a case, change Nts, Nt_ts and/or Nfm.
        end
        % Time-expansion error estimation.
        % Computing the extended source term at the test time-points:
        if Nt_ts_is_even
            if there_is_ih
                % The extended "source term" at the midpoint of the time-step is just the source term: 
                s_ext(:, texp_er_i) = s(:, texp_er_i);
            end
            % If there's no source term, it remains 0 throughout the propagation.
        else
            s_ext(:, texp_er_i) = Gdiff_op(Ulast2(:, texp_er_i), t(texp_er_i), Ulast2(:, tmidi), t(tmidi), varargin{:});
            if there_is_ih
                s_ext(:, texp_er_i) = s_ext(:, texp_er_i) + s(:, texp_er_i);
            end
        end
        % Computing the additional divided differences for the test-points:
        Cnewton(:, texp_er_i) = new_divdif(t_2ts(1:total_Nt_ts)*4/Tts, s_ext(:, texp_er_i), diagonal);
        normUnew = norm(Unew(:, Nt_ts));
        normCnewton_t_er = norm(Cnewton(:, Nt_ts + 1));
        texp_error_cheap = Etexp_factor_cheap*normCnewton_t_er/normUnew;
        if Nt_ts_is_even
            texp_error_exact = Etexp_factor_even*normCnewton_t_er/normUnew;
        elseif options.texp_er_odd
            % The odd estimation is computed only if this is specified in the
            % options structure.
            texp_error_exact =...
                Etexp_factor_odd*norm(4*Cnewton(:, Nt_ts + 2) - ...
                                     Tts*Gop(Ulast(:, tmidi), t(tmidi), Cnewton(:, Nt_ts + 1), varargin{:}))/normUnew;
        end
        if options.display_mode && tol_mode
            if conv_error>tol_ts
                fprintf('Warning: The estimated error of the iterative process (%d) is larger than the requested tolerance.\nThe solution might be inaccurate (in time step No. %d).\n', conv_error, tsi)
            end
            if Nt_ts_is_even || options.texp_er_odd
                if texp_error_exact>tol_ts 
                    fprintf('Warning: The estimated error of the time expansion (%d) is larger than the requested tolerance.\nThe solution might be inaccurate (in time step No. %d).\n', texp_error_exact, tsi)
                end
            elseif texp_error_cheap>tol_ts 
                fprintf('Warning: The estimated error of the time expansion (%d) is larger than the requested tolerance.\nThe solution might be inaccurate (in time step No. %d).\n', texp_error_cheap, tsi)
            end
            if fUerror>tol_ts
                fprintf('Warning: The estimation of the error resulting from the function of matrix (%d) is larger than the requested tolerance.\nThe solution might be inaccurate (in time step No. %d).\n', fUerror, tsi)
            end
        end        
        est_errors.texp_cheap = est_errors.texp_cheap + texp_error_cheap;
        if Nt_ts_is_even || options.texp_er_odd
            est_errors.texp_exact = est_errors.texp_exact + texp_error_exact;
        end
        est_errors.fm = est_errors.fm + fUerror;
        est_errors.conv = est_errors.conv + conv_error;
        est_errors.conv_texp = est_errors.conv_texp + conv_error_texp;
        est_errors.conv_fm = est_errors.conv_fm + conv_error_fm;
        if options.conv_er_cheb
            est_errors.conv_cheb = est_errors.conv_cheb + conv_error_cheb;
        end
        if nargout>4
            history.texp_error_cheap(tsi) = texp_error_cheap;
            if Nt_ts_is_even || options.texp_er_odd
                history.texp_error_exact(tsi) = texp_error_exact;
            end
            history.fm_error(tsi) = fUerror;
            history.reldif(tsi) = reldif;
            history.conv_error(tsi) = conv_error;
            history.conv_error_texp(tsi) = conv_error_texp;
            history.conv_error_fm(tsi) = conv_error_fm;
            if options.conv_er_cheb
                history.conv_error_cheb(tsi) = conv_error_cheb;
            end
            history.niter(tsi) = niter;
        end
        if tsi ~= 1
            allniter = allniter + niter;
        else
            % The iterations of the first time-step aren't summed in
            % allniter.
            if Arnoldi
                matvecs.propagation = niter*(Nt_ts + Ns_ext_i*Gdiff_matvecs + Nfm);
            else
                matvecs.propagation = niter*(Nt_ts + Ns_ext_i*Gdiff_matvecs + Nfm - 1);
            end
        end
        % Computation of the solution at the tgrid points.
        % Finding the indices of the tgrid points within the time step (the indices of the points
        % to be computed are between tgrid_lowi and tgrid_upi:
        while tgrid_upi<Nt && (t(Nt_ts) - tgrid(tgrid_upi + 1))*direction>abs(t(Nt_ts))*eps*10
            tgrid_upi = tgrid_upi + 1;
        end
        % Calculating the solution at the tgrid points: 
        if tgrid_lowi<=tgrid_upi
            timeMout = maketM(tgrid(tgrid_lowi:tgrid_upi) - t(1), Nt_ts);
            if Arnoldi
                U(:, tgrid_lowi:tgrid_upi) = Ufrom_vArnoldi(v_vecs, timeMout, Upsilon, RvKr, samplingp, capacity);
            else
                Ccheb_f_out = f_chebC(tgrid(tgrid_lowi:tgrid_upi) - t(1));
                U(:, tgrid_lowi:tgrid_upi) = Ufrom_vCheb(v_vecs, timeMout, Vcheb, Ccheb_f_out);
            end
            tgrid_lowi = tgrid_upi + 1;
        end
        % If one of the points in tgrid coincides with the point of the
        % propagation grid:
        if abs(t(Nt_ts) - tgrid(tgrid_upi + 1))<=abs(t(Nt_ts))*eps*10
            tgrid_upi = tgrid_upi + 1;
            U(:, tgrid_upi) = Unew(:, Nt_ts);
            tgrid_lowi = tgrid_upi + 1;
        end
        if nargout>4 && ~options.save_memory
            history.U(:, ((tsi - 1)*(Nt_ts - 1) + 2):(tsi*(Nt_ts - 1) + 1)) = Unew(:, 2:Nt_ts);
        end
        % The new guess is an extrapolation of the solution within the
        % previous time step:
        Uguess(:, 1) = Unew(:, Nt_ts);
        if Arnoldi
            Uguess(:, 2:total_Nt_ts) = Ufrom_vArnoldi(v_vecs, timeMnext, Upsilon, RvKr, samplingp, capacity);
        else
            Uguess(:, 2:total_Nt_ts) = Ufrom_vCheb(v_vecs, timeMnext, Vcheb, Ccheb_f_next);
        end
        % Divergence investigation:
%         fprintf('New guess norm:\n')
%         vecnorm(Uguess)
        if there_is_ih
            s(:, 1) = s(:, Nt_ts);
        end
    end
    if Nt_ts_is_even || options.texp_er_odd
        est_errors.total = est_errors.texp_exact + est_errors.fm + est_errors.conv;
    else
        est_errors.total = est_errors.texp_cheap + est_errors.fm + est_errors.conv;
    end
    if tol_mode
        if Nt_ts_is_even || options.texp_er_odd
            if est_errors.texp_exact>tol
                fprintf('\nWarning: The total estimated error of the time expansion (%d) is larger than the requested tolerance.\nThe solution might be inaccurate.\n', est_errors.texp_exact)
            end
        elseif est_errors.texp_cheap>tol
            fprintf('\nWarning: The total estimated error of the time expansion (%d) is larger than the requested tolerance.\nThe solution might be inaccurate.\n', est_errors.texp_cheap)
        end
        if est_errors.fm>tol
            fprintf('\nWarning: The total estimated error resulting from the function of matrix computation (%d) is larger than the requested tolerance.\nThe solution might be inaccurate.\n', est_errors.fm)
        end
        if est_errors.conv>tol
            fprintf('\nWarning: The total estimated error resulting from the iterative process (%d) is larger than the requested tolerance.\nThe solution might be inaccurate.\n', est_errors.conv)
        end
        if est_errors.total>tol
            fprintf('\nWarning: The total estimated error (%d) is larger than the requested tolerance.\nThe solution might be inaccurate.\n', est_errors.total)
        end
    end
    if nargout>1
        mniter = allniter/(Nts - 1);
    end
    if nargout>2
        % The number of Hamiltonian operations required for the propagation:
        if Arnoldi
            matvecs.propagation = matvecs.propagation + allniter*(Nt_ts + Ns_ext_i*Gdiff_matvecs + Nfm);
        else
            matvecs.propagation = matvecs.propagation + allniter*(Nt_ts + Ns_ext_i*Gdiff_matvecs + Nfm - 1);
        end
        % The number of Hamiltonian operations required for the error estimation:
        if options.conv_er_cheb
            matvecs.error = Nts*Gdiff_matvecs*Ns_ext_i;
            % The situation that this is more than 0 isn't recommended.
        else
            matvecs.error = Nts*Gdiff_matvecs;
        end
        if Nt_ts_is_even || ~options.texp_er_odd
            matvecs.error = matvecs.error + Nts*Gdiff_matvecs;
        else
            matvecs.error = matvecs.error + Nts*(2*Gdiff_matvecs + 1);
        end
        matvecs.total = matvecs.propagation + matvecs.error;
    end
    if nargout>4
        if Nt_ts_is_even || options.texp_er_odd
            history.total_error = history.texp_error_exact + history.fm_error + history.conv_error;
        else
            history.total_error = history.texp_error_cheap + history.fm_error + history.conv_error;
        end
    end
    
    %%% Nested functions: %%%
    
    function result = f_fun(z, t)
    % The function f(z, t):
        Ntp = length(t);
        Nz = length(z);
        zt = z*t;
        % Condition for estimating if f(z, t) should be computed directly or by
        % a "tail" of a Taylor expansion (see supplementary material):
        is_big = factorialNt_ts*eps./abs(zt.^(Nt_ts)) < options.tol_f;
        %is_big = factorialNt_ts./abs(zt.^(Nt_ts)) < 1;
        result = ones(Nz, Ntp);
        % First, we compute f(z, t)/(t^Nt_ts), which is a function of zt.
        % A direct computation for large arguments:
        result(is_big) = exp(zt(is_big));
        for polynomi = 1:Nt_ts
            result(is_big) = polynomi*(result(is_big) - 1)./zt(is_big);
        end
        % Computation by a Taylor form for small arguments:
        is_not_converged = ~is_big;
        term = double(is_not_converged);
        polydeg = 1;
        while max(max(is_not_converged))
            term(is_not_converged) = zt(is_not_converged).*term(is_not_converged)/(polydeg + Nt_ts);
            result(is_not_converged) = result(is_not_converged) + term(is_not_converged);
            polydeg = polydeg + 1;
            is_not_converged(is_not_converged) = abs(term(is_not_converged))./abs(result(is_not_converged)) > eps;
        end
        % Obtaining the required function f(z,t):
        result = result.*((ones(Nz, 1)*t).^Nt_ts);
    end

    function Rv = getRv(v_vecs, Hessenberg, samplingp, capacity)
        % Obtaining the R_n(G)v_vecs(:, Nt_ts + 1) represented in the
        % Krylov space, where the R_n(z)'s are the Newton basis polynomials,
        % with samplingp as the sampling points.
        Rv = zeros(Nfm + 1, Nfm + 1);
        % Rv(:, 1) is v_vecs(:, Nt_ts + 1) in the Krylov space.
        Rv(1, 1) = norm(v_vecs(:, Nt_ts + 1));
        for spi = 1:Nfm
            % Rv(:, spi) belongs to a Krylov space of dimension spi, and
            % the terms with indices larger than spi vanish.
            Rv(1:(spi + 1), spi + 1) = (Hessenberg(1:(spi + 1), 1:spi)*Rv(1:spi, spi) - samplingp(spi)*Rv(1:(spi + 1), spi))/capacity;
        end
    end
    
    function [U, f_abs_error, fUerror] = Ufrom_vArnoldi(v_vecs, timeM, Upsilon, RvKr, samplingp, capacity)
        % The function computes the solution from the v vectors for the
        % Arnoldi algorithm at all time points specified by the transpose
        % Vandermonde matrix timeM.
        % Input:
        % v_vecs: The v vectors
        % timeM: The matrix of the t powers for the required time-points
        % Upsilon: The Krylov space orthonormalized vectors
        % RvKr: The vectors computed by the function getRv
        % samplingp: The sampling points for the Newton expansion in the
        % Krylov space
        % capacity: The capacity of the approximation domain.
        % Output:
        % U: The solution in the required time points
        % f_abs_error: The estimated absolute error of the computation of f_fun(G,t(Nt_ts))v_vecs(:, Nt_ts + 1)
        % fUerror: The estimated relative error of U resulting from the
        % computational error of f(G,t(Nt_ts))v_vecs(:, Nt_ts + 1)
        
        % The required time-points are given by their first power:
        tp = timeM(2, :);
        % Computing the divided differences for the Newton expansion of f_fun, for
        % all time-points:
        fdvd_ts = divdif(samplingp.'/capacity, f_fun(samplingp, tp).').';
        % The computation of the Newton expansion of f_fun(G, tp)v_vecs(:, Nt_ts + 1)
        % in the Krylov space, for all tp:
        fGv_kr = RvKr(1:Nfm, :)*fdvd_ts;
        % Computing the solution in all tp: 
        U = v_vecs(:, 1:(end-1))*timeM + Upsilon(:, 1:Nfm)*fGv_kr;        
        if nargout>1
            % The absolute error:
            f_abs_error = abs(fdvd_ts(Nfm + 1, Nt_ts - 1))*norm(RvKr(:, Nfm + 1));
%             % The relative error:
%             f_error = f_abs_error/norm(fGv_kr(:, Nt_ts - 1));
            % The relative error of U, resulting from this computation:
            fUerror = f_abs_error/norm(U(:, Nt_ts - 1));
        end
    end

    function [U, fUerror] = Ufrom_vCheb(v_vecs, timeM, Vcheb, Ccheb_f)
        % The function computes the solution from the v vectors for the
        % Arnoldi algorithm at all time points specified by the transpose
        % Vandermonde matrix timeM.
        % Input:
        % v_vecs: The v vectors
        % timeM: The matrix of the t powers for the required time-points
        % Vcheb: The T_k(\tilde{G})v_{Nt_ts} vectors, k=0,1,...,Nt_ts-1, in
        % separate columns
        % Ccheb_f: The Chebyshev coefficients of \tilde{f}_{Nt_ts}(z, t) in the
        % required time-points specified by timeM, as computed by the
        % function f_chebC.
        % Output:
        % U: The solution in the required time points
        % fUerror: The estimated relative error of U resulting from the
        % computational error of f(G,t(Nt_ts))v_vecs(:, Nt_ts + 1)       
        fGv = Vcheb*Ccheb_f;
        U = v_vecs(:, 1:Nt_ts)*timeM + fGv;
        if nargout>1
            f_abs_error = f_scalar_error*norm(v_vecs(:, Nt_ts + 1));
            fUerror = f_abs_error/norm(U(:, Nt_ts - 1)); %%% Changed %%%
        end
    end

    function Ccheb_f = f_chebC(t)
    % The function computes the Chebyshev coefficients of f(z,t), where t serves as a
    % parameter. t represents a row vector of time-values.
        zsamp_cheb = cos(((1:Nfm).'*2 - 1)*pi/(2*Nfm));
        zsamp = 0.5*(zsamp_cheb*(max_ev - min_ev) + max_ev + min_ev);
        f_zt = f_fun(zsamp, t);
        Ccheb_f = chebcM(f_zt);
    end



end

%%%% Sub-functions: %%%

function timeM = maketM(t, Nt_ts)
% Computation of the matrix timeM of time Taylor polynomials.
% timeM(i, j) = t(j)^(i - 1)
    Nt = length(t);
    timeM = zeros(Nt_ts, Nt);
    timeM(1, :) = ones(1, Nt);
    for vi = 2:Nt_ts
        timeM(vi, :) = t.*timeM(vi - 1, :);
    end
end

function capacity = get_capacity(sp, testp)
    % Computation of the capacity of the eigenvalue domain:
    capacity = 1;
    sp_comp = sp(sp ~= testp);
    Nsp = length(sp_comp);
    for zi = 1:Nsp
        capacity = capacity*abs(testp - sp_comp(zi))^(1/Nsp);
    end
end


function Uguess = guess0(ui, Np)
    % The function returns the zeroth order approximation for the guess of the
    % first step:
    Uguess = ui*ones(1, Np);
end