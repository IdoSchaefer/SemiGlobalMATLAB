function lfuns = conv_lfuns(p_max, xivals)
% The function computes the values of l^{(p)}(xi) at values xivals, for all
% p up to the maximal iteration number p_max.
% Input:
% p_max: The maximal iteration number
% xivals: A row vector of the xi values to be computed
    C1st = convC1st(p_max);
    Mxi = xivals.^((0:(2*p_max + 1)).');
    lfuns = C1st*Mxi;
end