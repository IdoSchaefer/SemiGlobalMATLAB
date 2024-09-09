function C = convCfm(MplusK, p_max)
% The function computes the coefficients of (xi + 1)^j of h(xi) for the convergence
% process of the semi-global propagator for the function of matrix error.
% MplusK equals M + K. The function returns the coefficients of all iterations up to p_max.
% C: A matrix; contains the coefficients of the p'th iteration in the 
% (p + 1)'th row, where different columns represent the coefficients of
% different polynomial orders.
% Called function: recurrence_coefs.m
    Cinit = zeros(1, MplusK + 2*p_max + 1);
    Cinit(MplusK + 1) = 1;
    C = recurrence_coefs(Cinit, p_max, @f_rec);

    function newC = f_rec(previousC, p)
        newC = zeros(1, length(previousC));
        n = 0:(MplusK + 2*p);
        prev_i = 1:(MplusK + 2*p - 1);
        % The index of the order j coefficient is j+1.
        newC(1) = sum(previousC(prev_i).*(n(prev_i) + 4)./((n(prev_i) + 1).*(n(prev_i) + 2)))/2;
        newC(2) = -1.5*previousC(1);
        ni = 3:(MplusK + 2*p + 1);
        newC(ni) = (previousC(ni - 2) - 1.5*previousC(ni - 1))./n(ni);
    end

end
