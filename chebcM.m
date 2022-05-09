function C = chebcM(fM)
% The function computes the Chebychev coefficients of a
% function sampled at the Chebychev points.
% fM is a matrix that contains the sampled values of several functions in
% its columns.
% The Chebyshev coefficiets of each function are the corresponding columns
% of C.
    N = size(fM, 1);
    C = dct(fM)/sqrt(N);
    C(2:N, :) = C(2:N, :)*sqrt(2);
end