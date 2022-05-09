function c = chebc(f, leftb, rightb, N)
% The function returns the Chebychev coefficients of the function f in a given domain.
% f is a function handle. leftb and rightb are the boundaries of the
% domain. N is the number of Chebychev coefficients.
    % The Chebyshev points in the Chebyshev polynomial domain, [-1, 1]:
    xcheb = cos(((1:N)*2 - 1)*pi/(2*N));
    % The Chebyshev points transformed to the approximation domain,
    % [leftb, rightb]:
    x = 0.5*(xcheb*(rightb - leftb) + rightb + leftb);
    c = dct(f(x))/sqrt(N);
    c(2:N) = c(2:N)*sqrt(2);
end