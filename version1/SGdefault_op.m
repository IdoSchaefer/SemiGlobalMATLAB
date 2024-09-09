function options = SGdefault_op
% The function sets the default options for the SemiGlobal1.m code.
    % The maximal allowed number of iterations for all time-steps excluding
    % the first:
    options.Niter = 7;
    % The maximal allowed number of iterations for the first time-step:
    options.Niter1st = 10;
    % The field display_mode is a boolean; true means that all warnings are 
    % displayed, even during propagation, including warnings of error that
    % exceeds the tolerance. false means that only warnings of possibility
    % of failure of the entire process are displayed during propagation.
    % Warnings of error that exceeds the tolerance are displayed only at
    % the end of the propagation.
    options.display_mode = false;
    % The maximal Nt_ts:
    options.Nt_ts_max = 13;
    % The maximal Nfm:
    options.Nfm_max = 15;
    % The tolerance of the f_fun computation; empty means that it is set to
    % the default value: The tolerance per time-step if tol is specified,
    % and eps otherwise.
    options.tol_f = [];
    % The tolerance of the iterative process in the first time-step, if it 
    % is desired to treat it differently than the other time-steps; applies
    % also when tol isn't specified by the user, and the manual mode is
    % applied for the other time-steps.
    % empty means that it is set to the default value of the tolerance per
    % time-step (tol_ts).
    options.tol1st = [];
    % The field conv_er_cheb is a boolean. true means that the convergence
    % error estimation by integration at the Chebyshev time-points is performed.
    % false means the opposite. The default empty variable
    % means that this is determined by the value of Gdiff_matvecs.
    options.conv_er_cheb = [];
    % Indicates if the exact time-expansion error estimation for odd Nt_ts
    % is computed, or the cheap and inaccurate estimation only.
    options.texp_er_odd = true;
    % true means that the propagation history of the solution U isn't saved in
    % history.U.
    options.save_memory = false;
end