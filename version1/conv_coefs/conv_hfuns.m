function hfuns = conv_hfuns(MplusK, p_max, xivals)
% The function computes the values of h^{(p)}(xi) at values xivals, for all
% p up to the maximal iteration number p_max.
% Input:
% MplusK: The value of M+K
% p_max: The maximal iteration number
% xivals: A row vector of the xi values to be computed
    Cfm = convCfm(MplusK, p_max);
    Mxi_plus1 = (xivals + 1).^((0:(MplusK + 2*p_max)).');
    hfuns = Cfm*Mxi_plus1;
end