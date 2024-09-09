function [U, mniter, matvecs, max_errors, history] = SemiGlobal(Gop, Gdiff_op, Gdiff_matvecs,...
    ihfun, ev_domain, ui, tgrid, Nts, Nt_ts, Nfm, tol, Niter, Niter1st, display_mode, varargin)
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
% tgrid: the time grid of the desired output solution, U. Should be ordered in an
% increasing order for a forward propagation, and a decreasing order for a
% backward propagation.
% Nts: the number of time steps of the propagation. 
% Nt_ts: the number of interior Chebyshev time points in each time-step, used during
% the computational process (M in the article).
% Nfm: the number of expansion terms for the computation of the function of
% matrix (K in the article).
% tol: the desired tolerance of the convergence (epsilon in the article)
% Niter: The maximal allowed number of iterations for all time-steps
% excluding the first
% Niter1st: The maximal allowed number of iterations for the first
% time-step
% display_mode: A boolean variable. true means that warnings are displayed
% during the propagation. false means that warnings are displayed only
% before and after the propagation.
%
% Output:
% U: contains the solution at the time points specified in tgrid, in separate columns.
% mniter: the mean number of iteration for a time step, excluding the first step,
% which typically requires more iterations. Usually, should be 1 for ideal
% efficiency of the algorithm.
% matvecs: the number of Hamiltonian operations.
% max_errors: A structure which contains the maximal estimated errors: 
%   max_errors.texp: The maximal estimated error of U, resulting from the time-expansions in
%   each time-step
%   max_errors.fU: The maximal estimated error of U, resulting from the
%   computation of the function of matrix
%   max_errors.f = The maximal estimated error of the function of matrix
%   computation itself
%   max_errors.conv: The maximal estimated convergence error
% history: A structure which contains the history of propagation:
%   history.t: The time grid of propagation
%   history.U: The solution at history.t
%   history.texp_error: The estimations for the error of U, resulting from the 
%       time-expansion, for all time-steps
%   history.f_error: The estimations for the error of the computation of the
%       function of matrix for all time-steps (for the Arnoldi approximation only)
%   history.fUerror: The estimations for the error of U, resulting from the computation of the
%       function of matrix, for all time-steps
%   history.conv_error: The estimations for the convergence errors for all time-steps
%   history.niter: The number of iterations for each time-steps 

% Author: Ido Schaefer
% ido.schaefer@mail.huji.ac.il

%   Setting defaults:
    if nargin<14
        display_mode = true;
    end
    if nargin<12
        Niter = 10;
        Niter1st = 16;
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
    % The index of the middle term in the time step:
    tmidi = round(Nt_ts/2);
    Nu = length(ui);
    U = zeros(Nu, Nt);
    U(:, 1) = ui;
    % The Chebyshev points for expansion in time, in the domain in which the Chebyshev expansion
    % is defined: [-1 1]
    tcheb = -cos(((1:Nt_ts) - 1)*pi/(Nt_ts-1));
    % The Chebyshev points for expansion in time, in the domain of the time
    % variable:
    t_ts = 0.5*(tcheb + 1)*Tts;
    % Test point for error estimation:
    test_tpoint = (t_ts(tmidi) + t_ts(tmidi + 1))/2;
    % The interior time points of the current time step for interpolation,
    % and the next time step for extrapolation of the guess solution in the next step:
    t_2ts = [t_ts, test_tpoint, Tts + t_ts(2:Nt_ts), Tts + test_tpoint];    
    % The full propagation grid:
    propagation_grid = [kron((tinit:Tts:(tf - Tts)), ones(1, Nt_ts)) + kron(ones(1, Nts), [t_ts(1:(Nt_ts - 1)), test_tpoint]), tf];
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
        Ccheb_f_comp = f_chebC(t_2ts(2:(2*Nt_ts + 1)));
        Ccheb_f_ts = Ccheb_f_comp(:, 1:Nt_ts);
        Ccheb_f_next = Ccheb_f_comp(:, (Nt_ts + 1):(2*Nt_ts));   
        % Computing evenly spaced sampling points for the error test of the
        % function of matrix error:
        dz_test = (max_ev - min_ev)/(Nu + 1);
        ztest = min_ev + dz_test*(1:Nu).';
        fztest = f_fun(ztest, Tts);
        f_error = max(abs(chebc2result(Ccheb_f_ts(:, Nt_ts - 1), ev_domain, ztest) - fztest)./abs(fztest));
        if nargout>3
            max_errors.f = f_error;
        end
        if f_error>1e-5
            fprintf('Warning: The estimated error of the computation of the function of matrix (%d) is larger than 1e-5.\nInstability in the propagation process may occur.\n', f_error)
        end
    end
    if nargout>4
        history.U = zeros(Nu, Nts*(Nt_ts - 1) + 1);
        history.U(:, 1) = ui;
        history.t = [kron(tinit:Tts:(T - Tts), ones(1, Nt_ts - 1)), tf] + [kron(ones(1, Nts), t_ts(1:(Nt_ts - 1))), 0];
        history.texp_error = zeros(1, Nts);
        if Arnoldi
            % If the Arnoldi approximation is used, the error of the
            % function of matrix computation is different for each
            % time-step:
            history.f_error = zeros(1, Nts);
        end
        history.fUerror = zeros(1, Nts);
        history.conv_error = zeros(1, Nts);
        history.niter = zeros(1, Nts);
    end
    % Computing the matrix of the time Taylor polynomials.
    % timeMts contains the points in the current time step, and timeMnext
    % contains the points in the next time step:
    timeMcomp = maketM(t_2ts(2:(2*Nt_ts + 1)), Nt_ts);
    timeMts = timeMcomp(:, 1:Nt_ts);
    timeMnext = timeMcomp(:, (Nt_ts + 1):(2*Nt_ts));    
    % Computing the coefficients of the transformation from the Newton
    % interpolation polynomial terms, to a Taylor like form:
    Cr2t = r2Taylor4(t_ts, Tts);
    % The extended "inhomogeneous" vectors:
    s_ext = zeros(Nu, Nt_ts + 1);
    % The indices for the computation of s_ext (excluding tmidi):
    s_ext_i = [1:(tmidi - 1), (tmidi + 1):(Nt_ts + 1)];
    % The v vectors are defined recursively, and contain information about
    % the time dependence of the s_ext vectors:
    v_vecs = zeros(Nu, Nt_ts + 1);
    % If there is no inhomogeneous term in the equation, ihfun is empty, and there_is_ih = false.
    % If there is, there_is_ih = true.
    there_is_ih = ~isempty(ihfun);
    if there_is_ih
        s = zeros(Nu, Nt_ts + 1);
        s(:, 1) = ihfun(tinit, varargin{:});
    end
    % The 0'th order approximation is the first guess for the first time step.
    % Each column represents an interior time point in the time step:
    Uguess = guess0(ui, Nt_ts + 1);
    Unew = zeros(Nu, Nt_ts + 1);
    allniter = 0;
    % These variables are used to determine which points in tgrid are in the computed time-step.  
    tgrid_lowi = 2;
    tgrid_upi = 1;
    % The maximal errors:
    max_texp_error = 0;
    if Arnoldi        
        max_f_error = 0;
    end
    max_fUerror = 0;
    max_conv_error = 0;
    for tsi = 1:Nts
        % The time of the interior time points within the time-step:
        t = propagation_grid((tsi - 1)*Nt_ts + [(1:(Nt_ts - 1)), Nt_ts + 1, Nt_ts]);
        % The first guess for the iterative process, for the convergence of the u
        % values. Each column represents an interior time point in the time step:
        Ulast = Uguess;
        Unew(:, 1) = Ulast(:, 1);
        v_vecs(:, 1) = Ulast(:, 1);
        if there_is_ih
            % Computing the inhomogeneous term:
            s(:, 2:(Nt_ts + 1)) = ihfun(t(2:(Nt_ts + 1)), varargin{:});
        end
        % Starting an iterative process, until convergence:
        niter = 0;
        reldif = tol + 1;
        while (reldif>tol && ((tsi>1 && niter<Niter) || (tsi == 1 && niter<Niter1st)))
            % Calculation of the inhomogeneous s_ext vectors. Note that
            % s_ext(:, tmidi) is equivalent to s(:, tmidi), and therefore not
            % calculated.
            s_ext(:, s_ext_i) = Gdiff_op(Ulast(:, s_ext_i), t(s_ext_i), Ulast(:, tmidi), t(tmidi), varargin{:});
            % If there is an inhomogeneous term, we add it to the s_ext
            % vectors:
            if there_is_ih
                s_ext(:, tmidi) = s(:, tmidi);
                s_ext(:, s_ext_i) = s_ext(:, s_ext_i) + s(:, s_ext_i);
            end
            % Calculation of the coefficients of the form of Taylor
            % expansion, from the coefficients of the Newton
            % interpolation at the points t_ts.
            % The divided differences are computed by the function divdif.
            % For numerical stability, we have to transform the time points
            % in the time step, to points in an interval of length 4:
            Cnewton = divdif(t_ts*4/Tts, s_ext(:, 1:Nt_ts));
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
                if nargout>3
                    max_errors.texp = max_texp_error;
                    max_errors.fU = max_fUerror;
                    if Arnoldi
                        max_errors.f = max_f_error;
                    end
                    max_errors.conv = max_conv_error;
                end
                return
            end
            if Arnoldi
                % Creating the Krylov space by the Arnodi iteration procedure,
                % in order to approximate f(G,t)v_vecs(:, Nt_ts + 1):
                [Upsilon, Hessenberg] = createKrop(@(v) Gop(Ulast(:, tmidi), t(tmidi), v, varargin{:}), v_vecs(:, Nt_ts + 1), Nfm);
                % Obtaining eigenvalues of the Hessenberg matrix:
                eigval = eig(Hessenberg(1:Nfm, 1:Nfm));
                % The test point is the average point of the eigenvalues:
                avgp = sum(eigval)/Nfm;
                samplingp = [eigval; avgp];
                capacity = get_capacity(eigval, avgp);
                % Obtaining the expansion vectors for a Newton approximation of f(G,t)v_vecs(:, Nt_ts + 1)
                % in the reduced Krylov space:
                RvKr = getRv(v_vecs, Hessenberg, samplingp, capacity);
                % Calculation of the solution in all time points
                % within the time step:
                [Unew(:, 2:(Nt_ts + 1)), f_error, fUerror] = Ufrom_vArnoldi(v_vecs, timeMts, Upsilon, RvKr, samplingp, capacity);
            else
                % Employing a Chebyshev approximation for the function of
                % matrix computation.
                % Vcheb is a matrix. Contains the vectors:
                % Tn(G(Ulast(:,tmidi), t(tmidi)))*v_vecs(: ,Nt_ts + 1),  n = 0, 1, ..., Nfm-1
                % where the Tn(z) are the Chebyshev polynomials.
                % The n'th vector is the (n+1)'th column of Vcheb.
                Vcheb = vchebMop(@(v) Gop(Ulast(:, tmidi), t(tmidi), v, varargin{:}), v_vecs(:, Nt_ts + 1), min_ev, max_ev, Nfm);
                % Calculation of the solution in all the time points
                % within the time step:
                [Unew(:, 2:(Nt_ts + 1)), fUerror] = Ufrom_vCheb(v_vecs, timeMts, Vcheb, Ccheb_f_ts);
            end
            % Error estimation for the time expansion.
            % A Newton interpolation of s_ext in the test time-point:
            s_ext_testpoint = Cnewton(:, Nt_ts);
            for ti = (Nt_ts-1):-1:1
                s_ext_testpoint = Cnewton(:, ti) + (4/Tts)*(test_tpoint - t_ts(ti))*s_ext_testpoint;
            end
            % The error estimation of the time-expansion:
            texpansion_error = norm(s_ext_testpoint - s_ext(:, Nt_ts + 1))*abs(Tts)/norm(Unew(:, Nt_ts + 1));
            % Convergence check of u:
            reldif = norm(Unew(:, Nt_ts) - Ulast(:, Nt_ts))/norm(Ulast(:, Nt_ts));
            Ulast = Unew;
            niter = niter + 1;
        end
        if display_mode && reldif>tol
            fprintf('Warning: The program has failed to achieve the desired tolerance in the iterative process (in time step No. %d).\n', tsi)
            % In such a case, change Nts, Nt_ts and/or Nfm.
        end
        if display_mode && texpansion_error>tol 
            fprintf('Warning: The estimated error of the time expansion (%d) is larger than the requested tolerance.\nThe solution might be inaccurate (in time step No. %d).\n', texpansion_error, tsi)
        end
        if display_mode && fUerror>tol 
            fprintf('Warning: The estimation of the error resulting from the function of matrix (%d) is larger than the requested tolerance.\nThe solution might be inaccurate (in time step No. %d).\n', fUerror, tsi)
        end
        if display_mode && Arnoldi && f_error>1e-5
            fprintf('Warning: The estimated error of the function of matrix (%d) is larger than 1e-5.\nInstability in the propagation process may occur (in time step No. %d).\n', f_error, tsi)
        end
        % Updating the maximal errors:
        if max_texp_error<texpansion_error
            max_texp_error = texpansion_error;
        end
        if Arnoldi && max_f_error<f_error
            max_f_error = f_error;
        end
        if max_fUerror<fUerror
            max_fUerror = fUerror;
        end
        if max_conv_error<reldif
            max_conv_error = reldif;
        end
        if nargout>4
            history.texp_error(tsi) = texpansion_error;
            if Arnoldi
                history.f_error(tsi) = f_error;
            end
            history.fUerror(tsi) = fUerror;
            history.conv_error(tsi) = reldif;
            history.niter(tsi) = niter;
        end
        if tsi ~= 1
            allniter = allniter + niter;
        else
            if Arnoldi
                matvecs = niter*(Nt_ts*(1 + Gdiff_matvecs) + Nfm);
            else
                matvecs = niter*(Nt_ts*(1 + Gdiff_matvecs) + Nfm - 1);
            end
        end
        % Computation of the solution at the tgrid points.
        % Finding the indices of the tgrid points within the time step (the indices of the points
        % to be computed are between tgrid_lowi and tgrid_upi:
        while tgrid_upi<Nt && (t(Nt_ts) - tgrid(tgrid_upi + 1))*direction>abs(t(Nt_ts))*eps*10
            tgrid_upi = tgrid_upi + 1;
        end
        % Calculating the solution in the tgrid points: 
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
        if nargout>4
            history.U(:, ((tsi - 1)*(Nt_ts - 1) + 2):(tsi*(Nt_ts - 1) + 1)) = Unew(:, 2:Nt_ts);
        end
        % The new guess is an extrapolation of the solution within the
        % previous time step:
        Uguess(:, 1) = Unew(:, Nt_ts);
        if Arnoldi
            Uguess(:, 2:(Nt_ts + 1)) = Ufrom_vArnoldi(v_vecs, timeMnext, Upsilon, RvKr, samplingp, capacity);
        else
            Uguess(:, 2:(Nt_ts + 1)) = Ufrom_vCheb(v_vecs, timeMnext, Vcheb, Ccheb_f_next);
        end
        if there_is_ih
            s(:, 1) = s(:, Nt_ts);
        end
    end
    if nargout>3
        max_errors.texp = max_texp_error;
        max_errors.fU = max_fUerror;
        if Arnoldi
            max_errors.f = max_f_error;
        end
        max_errors.conv = max_conv_error;
    end
    if max_texp_error>tol
        fprintf('\nWarning: The maximal estimated error of the time expansion (%d) is larger than the requested tolerance.\nThe solution might be inaccurate.\n', max_texp_error)
    end
    if max_fUerror>tol
        fprintf('\nWarning: The maximal estimated error resulting from the function of matrix computation (%d) is larger than the requested tolerance.\nThe solution might be inaccurate.\n', max_fUerror)
    end
    if Arnoldi && max_f_error>1e-5
        fprintf('\nWarning: The maximal estimated error of the function of matrix (%d) is larger than 1e-5.\nInstability in the propagation process is possible.\n', max_f_error)
    end
    if max_conv_error>tol
        fprintf('\nWarning: The maximal estimated error resulting from the iterative process (%d) is larger than the requested tolerance.\nThe solution might be inaccurate.\n', max_conv_error)
    end
    if nargout>1
        mniter = allniter/(Nts - 1);
    end
    if nargout>2
        if Arnoldi
            matvecs = matvecs + allniter*(Nt_ts*(1 + Gdiff_matvecs) + Nfm);
        else
            matvecs = matvecs + allniter*(Nt_ts*(1 + Gdiff_matvecs) + Nfm - 1);
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
        is_big = factorialNt_ts*eps./abs(zt.^(Nt_ts)) < tol;
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
    
    function [U, f_error, fUerror] = Ufrom_vArnoldi(v_vecs, timeM, Upsilon, RvKr, samplingp, capacity)
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
        % f_error: The estimated relative error of the computation of f_fun(G,t(Nt_ts))v_vecs(:, Nkr + 1)
        % fUerror: The estimated relative error of U resulting from the
        % computational error of f(G,t(Nt_ts))v_vecs(:, Nkr + 1)
        
        % The required time-points are given by their first power:
        tp = timeM(2, :);
        % Computing the divided differences for the Newton expansion of f_fun, for
        % all time-points:
        fdvd_ts = divdif(samplingp.'/capacity, f_fun(samplingp, tp).').';
        % The computation of the Newton expansion of f_fun(G, tp)v_vecs(:, Nkr + 1)
        % in the Krylov space, for all tp:
        fGv_kr = RvKr(1:Nfm, :)*fdvd_ts;
        % Computing the solution in all tp: 
        U = v_vecs(:, 1:(end-1))*timeM + Upsilon(:, 1:Nfm)*fGv_kr;        
        if nargout>1
            % The absolute error:
            f_error_abs = abs(fdvd_ts(Nfm + 1, Nt_ts - 1))*norm(RvKr(:, Nfm + 1));
            % The relative error:
            f_error = f_error_abs/norm(fGv_kr(:, Nt_ts - 1));
            % The relative error of U, resulting from this computation:
            fUerror = f_error_abs/norm(U(:, Nt_ts - 1));
        end
    end

    function [U, fUerror] = Ufrom_vCheb(v_vecs, timeM, Vcheb, Ccheb_f)
        
        fGv = Vcheb*Ccheb_f;
        U = v_vecs(:, 1:(end-1))*timeM + fGv;
        if nargout>1
            fUerror = f_error*norm(fGv(:, end - 1))/norm(U(:, end - 1));
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