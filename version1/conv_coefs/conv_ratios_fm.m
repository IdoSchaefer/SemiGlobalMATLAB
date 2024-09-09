function all_eta = conv_ratios_fm(MplusKmax, p_max)
% The function computes the convergence ratios of the function of matrix
% extrapolation error. The ratios are computed for all 2=<M+K<=MplusKmax and for
% all iteration numbers up to p_max.
% The coefficients are defined for M>=2, so M+K = 1 isn't included in
% all_eta.
% all_eta: A matrix containing the covergence ratios; different M's
% are represented by different columns, and different iteration numbers by
% different rows.
    all_hfuns = zeros(p_max + 1, MplusKmax - 1);
    for MplusKi = 1:(MplusKmax - 1)
        all_hfuns(:, MplusKi) = conv_hfuns(MplusKi + 1, p_max, 1);
    end
    all_eta = all_hfuns(2:(p_max + 1), :)./all_hfuns(1:p_max, :);
end