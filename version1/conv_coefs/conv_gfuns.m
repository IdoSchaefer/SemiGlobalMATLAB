function gfuns = conv_gfuns(M, p_max, xivals)
% The function computes the values of g^{(p)}(xi) at values xivals, for all
% p up to the maximal iteration number p_max.
% Input:
% M: The M value
% p_max: The maximal iteration number
% xivals: A row vector of the xi values to be computed
    Ctexp = convCtexp(M, p_max);
    Mcheb = cheb_pols(2*xivals + 1, M + 2*p_max + 1);
    gfuns = Ctexp*Mcheb;
end