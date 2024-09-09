function data = SGdata(options)
    data.correct_texp = 2*conv_ratios_texp(options.Nt_ts_max, options.Niter);
    % For high Mmax, an alternative calculation of the convergence ratios
    % based on symbolic math can be included (conv_gfuns_sym.m), which is
    % more stable. This isn't urgent right now.
    data.correct_fm = 2*conv_ratios_fm(options.Nt_ts_max + options.Nfm_max, options.Niter);
    data.correct_1st = 2*conv_ratios_1st(options.Niter1st);
end