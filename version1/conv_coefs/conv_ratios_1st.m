function lambda = conv_ratios_1st(p_max)
% The function computes the convergence ratios of the first time-step.
% The ratios are computed for all iteration numbers up to p_max.
% lambda: A column vector containing the covergence ratios for all
% iteration numbers
    lfuns = conv_lfuns(p_max, 1);
    lambda = lfuns(2:(p_max + 1))./lfuns(1:p_max);
end