function all_gamma = conv_ratios_texp(Mmax, p_max)
% The function computes the convergence ratios of the time-expansion
% extrapolation error. The ratios are computed for all 2=<M<=Mmax and for
% all iteration numbers up to p_max.
% The coefficients are defined for M>=2, so M = 1 isn't included in
% all_gamma.
% all_gamma: A matrix containing the covergence ratios; different M's
% are represented by different columns, and different iteration numbers by
% different rows.
    all_gfuns = zeros(p_max + 1, Mmax - 1);
    for Mi = 1:(Mmax - 1)
        all_gfuns(:, Mi) = conv_gfuns(Mi + 1, p_max, 1);
    end
    all_gamma = all_gfuns(2:(p_max + 1), :)./all_gfuns(1:p_max, :);
end