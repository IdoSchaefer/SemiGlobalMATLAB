function [U, mniter, matvecs, max_errors, history] = SemiGlobalArnoldi_xp(K, V, Vtfun, ihfun, ui, x, tgrid, Nts, Nt_ts, Nkr, tol, Niter, Niter1st, display_mode, varargin)
% The program solves time-dependent Schroedinger equation for a
% time-dependent, nonlinear Hamiltonian, using the semi-global propagator 
% by Hillel Tal-Ezer.
% There is an option to include an inhomogeneos term.
% The program should be used if H0 is composed only from x diagonal and p
% diagonal elements, and the time dependent perturbation Vt is composed
% only from x diagonal elements.
%
% Input:
% K: the p diagonal element of the unperturbed Hamiltonian (kinetic
% energy).
% V: the x diagonal element of the unperturbed Hamiltonian (potential
% energy).
% Vtfun: A function handle of the form: @(u, x, t, more arguments). Returns
% the x vector of the time dependent, nonlinear perturbation: Vt(u, x, t).
% ihfun: A function handle of the form: @(t, more arguments). Returns the
% inhomogeneous term vector. If there is no inhomogeneous term, put []
% instead.
% Note: "more arguments" must be the same for the 2 functions. The arguments
% will be written in the place of "varargin" in this program.
% ui: the initial state vector.
% x: the x grid
% tgrid: the time grid of the desired output solution, U. Should be ordered in an
% increasing order for a forward propagation, and a decreasing order for a
% backward propagation.
% Nts: the number of time steps of the propagation. 
% Nt_ts: the number of interior Chebyshev time points in each time-step, used during
% the computational process (M in the article).
% Nkr: the dimension of the Krylov space in the Arnoldi procedure (K in the
% article).
% tol: the desired tolerance of the convergence (epsilon in the article)
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
%       function of matrix for all time-steps
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
    % In order to detect if the propagation is a forward or a backward
    % propagation:
    direction = sign(tgrid(2) - tgrid(1));
    Nt = length(tgrid);
    tinit = tgrid(1);
    tf = tgrid(Nt);
    % The length of the time interval of the whole propagation(can be negative):
    T = tf - tinit;
    % The length of the time step interval:
    Tts = T/Nts;
    % The index of the middle term in the time step:
    tmidi = round(Nt_ts/2);
    Nx = length(ui);
    U = zeros(Nx, Nt);
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
    if nargout>4
        history.U = zeros(Nx, Nts*(Nt_ts - 1) + 1);
        history.U(:, 1) = ui;
        history.t = [kron(0:Tts:(T - Tts), ones(1, Nt_ts - 1)), T] + [kron(ones(1, Nts), t_ts(1:(Nt_ts - 1))), 0];
        history.texp_error = zeros(1, Nts);
        history.f_error = zeros(1, Nts);
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
    % The extended "inhomogeneos" vectors:
    s_ext = zeros(Nx, Nt_ts + 1);
    % The indices for the computation od s_ext (excluding tmidi):
    s_ext_i = [1:(tmidi - 1), (tmidi + 1):(Nt_ts + 1)];
    % The v vectors are defined in recursively, and contain information about
    % the time dependence of the s_ext vectors:
    v_vecs = zeros(Nx, Nt_ts + 1);
    % If there is no inhomogeneous term in the equation, ihfun is empty, and there_is_ih = false.
    % If there is, there_is_ih = true.
    there_is_ih = ~isempty(ihfun);
    if there_is_ih
        s = zeros(Nx, Nt_ts + 1);
        s(:, 1) = ihfun(tinit, varargin{:});
    end
    % The 0'th order approximation is the first guess for the first time step.
    % Each column represents an interior time point in the time step:
    Uguess = guess_ts1(ui, Nt_ts + 1);
    allniter = 0;
    % These variables are used to determine which points in tgrid are in the computed time-step.  
    tgrid_lowi = 2;
    tgrid_upi = 1;
    % The maximal errors:
    max_texp_error = 0;
    max_f_error = 0;
    max_fUerror = 0;
    max_conv_error = 0;
    for tsi = 1:Nts
        % The time of the interior time points within the time-step:
        t = tinit + Tts*(tsi - 1) + [t_ts, test_tpoint];
        % The first guess for the iterative process, for the convergence of the u
        % values. Each column represents an interior time point in the time step:
        Ulast = Uguess;
        Unew(:, 1) = Ulast(:, 1);
        v_vecs(:, 1) = Ulast(:, 1);
        if there_is_ih
            % Computing the inhomogeneous term:
            for ti = 2:(Nt_ts + 1)
                s(:, ti) = ihfun(t(ti), varargin{:});
            end
        end
        % Starting an iterative process, until convergence:
        niter = 0;
        reldif = tol + 1;
        while (reldif>tol && ((tsi>1 && niter<Niter) || (tsi == 1 && niter<Niter1st)))
            % Vthalf is the x vector representing Vt(t = t(tmidi)): 
            Vthalf = Vtfun(Ulast(:, tmidi), x, t(tmidi), varargin{:});
            % Vts is the x vector representing the total potential energy 
            % operator at t = t(tmidi):
            Vts = V + Vthalf;
            % Calculation of the inhomogeneous s_ext vectors. Note that
            % s_ext(:, tmidi) is equivalent to s(:, tmidi), and therefore not
            % calculated.
            for ti = s_ext_i
                s_ext(:, ti) = -1i*(Vtfun(Ulast(:, ti), x, t(ti), varargin{:}) - Vthalf).*Ulast(:, ti);
            end
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
            Cnewton = divdif([t_ts, test_tpoint]*4/Tts, s_ext);
            % Calculating the Taylor like coefficients:
            Ctaylor = zeros(Nx, Nt_ts);
            Ctaylor(:, 1) = Cnewton(:, 1);
            for Newtoni = 2:Nt_ts
                Ctaylor(:, 1:Newtoni) = Ctaylor(:, 1:Newtoni) + Cnewton(:, Newtoni)*Cr2t(Newtoni, 1:Newtoni);
            end
            % Calculation of the v vectors:
            for polyi = 2:(Nt_ts+1)
                v_vecs(:, polyi) = (-1i*Hpsi(K, Vts ,v_vecs(:, polyi-1)) + Ctaylor(:, polyi-1))/(polyi - 1);
            end
            if ~min(isfinite(v_vecs(:, Nt_ts + 1)))
                % It means that the solution diverges.
                % In such a case, change Nts, Nt_ts and/or Ncheb.
                fprintf('\nError: The algorithm diverges (in time step No. %d).\n', tsi);
                mniter = allniter/(tsi - 1);
                if nargout>3
                    max_errors.texp = max_texp_error;
                    max_errors.fU = max_fUerror;
                    max_errors.f = max_f_error;
                    max_errors.conv = max_conv_error;
                end
                return
            end
            % Creating the Krylov space by the Arnodi iteration procedure,
            % in order to approximate f(G,t)v_vecs(:, Nt_ts + 1):
            [Upsilon, Hessenberg] = createKrop(@(u) Hpsi(K, Vts, u), v_vecs(:, Nt_ts + 1), Nkr);
%            [Upsilon, Hessenberg] = createKrop(@(u) -1i*Hpsi(K, Vts, u), v_vecs(:, Nt_ts + 1), Nkr);
            % Obtaining eigenvalues of the Hessenberg matrix:
            eigval = eig(Hessenberg(1:Nkr, 1:Nkr));
            % The test point is the average point of the eigenvalues:
            avgp = sum(eigval)/Nkr;
            samplingp = [eigval; avgp];
            capacity = get_capacity(eigval, avgp);
            % Obtaining the expansion vectors for a Newton approximation of f(-1i*H,t)v_vecs(:, Nt_ts + 1)
            % in the reduced Krylov space:
            RvKr = getRv;
            % Calculation of the solution in all time points
            % within the time step:
            [Unew(:, 2:(Nt_ts + 1)), f_error, fUerror] = Ufrom_v(timeMts);
            % Error estimation for the time expansion.
            % A Newton interpolation of s_ext in the test time-point:
            s_ext_testpoint = Cnewton(:, Nt_ts);
            for ti = (Nt_ts-1):-1:1
                s_ext_testpoint = Cnewton(:, ti) + (4/Tts)*(test_tpoint - t_ts(ti))*s_ext_testpoint;
            end
            % The error estimation of the time-expansion:
            texpansion_error = norm(s_ext_testpoint - s_ext(:, Nt_ts + 1))*abs(Tts)/norm(Unew(:, Nt_ts + 1));
            if display_mode && texpansion_error>tol 
                fprintf('Warning: The estimated error of the time expansion (%d) is larger than the requested tolerance.\nThe solution might be inaccurate (in time step No. %d).\n', texpansion_error, tsi)
            end
            if display_mode && fUerror>tol 
                fprintf('Warning: The estimation of the error resulting from the function of matrix (%d) is larger than the requested tolerance.\nThe solution might be inaccurate (in time step No. %d).\n', fUerror, tsi)
            end
            if display_mode && f_error>1e-5
                fprintf('Warning: The estimated error of the function of matrix (%d) is larger than 1e-5.\nInstability in the propagation process may occur (in time step No. %d).\n', f_error, tsi)
            end
            % Convergence check of u:
            reldif = norm(Unew(:, Nt_ts) - Ulast(:, Nt_ts))/norm(Ulast(:, Nt_ts));
            Ulast = Unew;
            niter = niter + 1;
        end
        if display_mode && niter == Niter
            display('Warning: The program has failed to achieve the desired tolerance in the iterative process.')
            % In such a case, change Nts, Nt_ts and/or Ncheb.
        end
        % Updating the maximal errors:
        if max_texp_error<texpansion_error
            max_texp_error = texpansion_error;
        end
        if max_f_error<f_error
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
            history.f_error(tsi) = f_error;
            history.fUerror(tsi) = fUerror;
            history.conv_error(tsi) = reldif;
            history.niter(tsi) = niter;
        end
        if tsi ~= 1
            allniter = allniter + niter;
        else
            matvecs = niter*(Nt_ts + Nkr);
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
            U(:, tgrid_lowi:tgrid_upi) = Ufrom_v(timeMout);
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
        Uguess(:, 2:(Nt_ts + 1)) = Ufrom_v(timeMnext);
        if there_is_ih
            s(:, 1) = s(:, Nt_ts);
        end
    end
    if nargout>3
        max_errors.texp = max_texp_error;
        max_errors.fU = max_fUerror;
        max_errors.f = max_f_error;
        max_errors.conv = max_conv_error;
    end
    if max_texp_error>tol
        fprintf('\nWarning: The maximal estimated error of the time expansion (%d) is larger than the requested tolerance.\nThe solution might be inaccurate.\n', max_texp_error)
    end
    if max_fUerror>tol
        fprintf('\nWarning: The maximal estimated error resulting from the function of matrix computation (%d) is larger than the requested tolerance.\nThe solution might be inaccurate.\n', max_fUerror)
    end
    if max_f_error>1e-5
        fprintf('\nWarning: The maximal estimated error of the function of matrix (%d) is larger than 1e-5.\nInstability in the propagation process is possible.\n', max_f_error)
    end
    if max_conv_error>tol
        fprintf('\nWarning: The maximal estimated error resulting from the iterative process (%d) is larger than the requested tolerance.\nThe solution might be inaccurate.\n', max_conv_error)
    end
    mniter = allniter/(Nts - 1);
    matvecs = matvecs + allniter*(Nt_ts + Nkr);
    
    %%% Nested functions: %%%
    
    function Rv = getRv
        % Obtaining the R_n(G)v_vecs(:, Nt_ts + 1) represented in the
        % Krylov space, where the R_n(z)'s are the Newton basis polynomials,
        % with samplingp as the sampling points.
        Rv = zeros(Nkr + 1, Nkr + 1);
        % Rv(:, 1) is v_vecs(:, Nt_ts + 1) in the Krylov space.
        Rv(1, 1) = norm(v_vecs(:, Nt_ts + 1));
        for spi = 1:Nkr
            % Rv(:, spi) belongs to a Krylov space of dimension spi, and
            % the terms with indices larger than spi vanish.
            Rv(1:(spi + 1), spi + 1) = (Hessenberg(1:(spi + 1), 1:spi)*Rv(1:spi, spi) - samplingp(spi)*Rv(1:(spi + 1), spi))/capacity;
        end
    end
    
    function [U, f_error, fUerror] = Ufrom_v(timeM)
        % The function computes the solution from the v vectors.
        % Input:
        % timeM: The matrix of the t powers for the required time-points
        % Output:
        % U: The solution in the required time points
        % f_error: The estimated relative error of the computation of f(-1i*H,t(Nt_ts))v_vecs(:, Nkr + 1)
        % fUerror: The estimated relative error of U resulting from the
        % computational error of f(-1i*H,t(Nt_ts))v_vecs(:, Nkr + 1)
        
        % The required time-points are given by their first power:
        tp = timeM(2, :);
        % Computing the divided differences for the Newton expansion of f_fun, for
        % all time-points:
        fdvd_ts = divdif(samplingp.'/capacity, f_fun(samplingp, tp, Nt_ts, tol).').';
        % The computation of the Newton expansion of f(-1i*H, tp)v_vecs(:, Nkr + 1)
        % in the Krylov space, for all tp:
        fHv_kr = RvKr(1:Nkr, :)*fdvd_ts;
        % Computing the solution in all tp: 
        U = v_vecs(:, 1:(end-1))*timeM + Upsilon(:, 1:Nkr)*fHv_kr;        
        if nargout>1
            % The absolute error:
            f_error_abs = abs(fdvd_ts(Nkr + 1, Nt_ts - 1))*norm(RvKr(:, Nkr + 1));
            % The relative error:
            f_error = f_error_abs/norm(fHv_kr(:, Nt_ts - 1));
            % The relative error of U, resulting from this computation:
            fUerror = f_error_abs/norm(U(:, Nt_ts - 1));
        end
    end

end

%%%% Sub-functions: %%%

function result = f_fun(z, t, Nt_ts, tol)
% The function f(-1i*z, t):
    Nt = length(t);
    Nz = length(z);
    minus_izt = -1i*z*t;
%    minus_izt = z*t;
    % Condition for estimating if f_fun(z, t) should be computed directly or by
    % a "tail" of a Taylor expansion:
    is_big = factorial(Nt_ts)*eps./abs(minus_izt.^(Nt_ts)) < tol;
    result = ones(Nz, Nt);
    % First, we compute f(-1i*z, t)/(t^Nt_ts), which is a function of zt.
    % A direct computation for large arguments:
    result(is_big) = exp(minus_izt(is_big));
    for polyi = 1:Nt_ts
        result(is_big) = polyi*(result(is_big) - 1)./minus_izt(is_big);
    end
    % Computation by a Taylor form for small arguments:
    is_not_converged = ~is_big;
    term = double(is_not_converged);
    polydeg = 1;
    while max(max(is_not_converged))
        term(is_not_converged) = minus_izt(is_not_converged).*term(is_not_converged)/(polydeg + Nt_ts);
        result(is_not_converged) = result(is_not_converged) + term(is_not_converged);
        polydeg = polydeg + 1;
        is_not_converged(is_not_converged) = abs(term(is_not_converged))./abs(result(is_not_converged)) > eps;
    end
    % Obtaining the required function f(-1i*z, t):
    result = result.*((ones(Nz, 1)*t).^Nt_ts);
end

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


function Uguess = guess_ts1(ui, Np)
    % The function returns the zeroth order approximation for the guess of the
    % first step:
    Uguess = ui*ones(1, Np);
end