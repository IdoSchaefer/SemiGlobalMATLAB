function c = chebcbv(fv)
% The function computes the Chebychev coefficients of a
% function sampled at the Chebychev points that include the boundary of the
% domain.
% fv is the vector of the sampled values.
    N = length(fv) - 1;
    c = sqrt(2/N)*dctI(fv);
    c([1, N+1]) = 0.5*c([1, N+1]);
end