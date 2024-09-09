function C = convC1st(p_max)
% The function computes the power coefficients of l(xi) for the convergence
% process of the semi-global propagator for the first time-step. The
% function returns the coefficients of all iterations up to p_max.
% C: A matrix; contains the coefficients of the p'th iteration in the 
% (p + 1)'th row, where different columns represent the coefficients of
% different orders of polynomials (the n'th column represents the coefficient
% of xi^(n-1)).
% Called function: recurrence_coefs.m
    Cinit = zeros(1, 2*p_max + 2);
    Cinit(2) = 1;
    C = recurrence_coefs(Cinit, p_max, @f_rec);
end

function newC = f_rec(previousC, p)
    newC = zeros(1, length(previousC));
    n = (p + 1):(2*p + 1);
    % The index of the order j coefficient is j+1.
    newC(n + 1) = (previousC(n - 1) - previousC(n)/2)./n;
end