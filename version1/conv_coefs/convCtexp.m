function C = convCtexp(M, p_max)
% The function computes the Chebyshev coefficients of g(xi) for the convergence
% process of the semi-global propagator for the time expansion error.
% The function returns the coefficients of all iterations up to p_max.
% C: A matrix; contains the coefficients of the p'th iteration in the 
% (p + 1)'th row, where different columns represent the coefficients of
% different polynomial orders.
% Called function: recurrence_coefs.m
    Cinit = zeros(1, M + 2*p_max + 2);
    if M==3
        Cinit(1) = 3/8;
    elseif mod(M, 2) == 0
        Cinit(1) = 4/((M^2 - 1)*(M - 3));
    else
        Cinit(1) = -4/((M^2 - 1)*(M - 3));
    end
    Cinit(M) = -1/(M - 1);
    Cinit(M + 2) = 1/(2*(M + 1));
    if M~=3
        % This accounts also for the case of M = 2:
        Cinit(abs(M - 3) + 1) = Cinit(abs(M - 3) + 1) + 1/(2*(M - 3));
    end
    Cinit = Cinit/4^M;
    C = recurrence_coefs(Cinit, p_max, @f_rec);

    function newC = f_rec(previousC, p)
        Nmax = length(previousC) - 1;
        newC = zeros(1, Nmax + 1);
        N = M + 2*p + 1;
        n = 0:N;
        prev_i = 4:(N - 1);
        % The index of the order j coefficient is j+1.
        newC(1) = (7/4)*previousC(1) + previousC(2)/6 - (35/48)*previousC(3)...
            - sum(previousC(prev_i).*(n(prev_i).^2 - 7)./((n(prev_i).^2 - 4).*(n(prev_i).^2 - 1)));
        newC(2) = -2*previousC(1) + previousC(2)/4 + previousC(3) - previousC(4)/4;
        newC(3) = previousC(1)/4 - previousC(2)/2 + previousC(4)/2 - previousC(5)/8;
        ni = 4:(N - 3);
        newC(ni) = ((previousC(ni - 2) - previousC(ni + 2))/4 - previousC(ni - 1) + previousC(ni + 1))./n(ni);
        newC(N - 2) = (previousC(N - 4)/4 - previousC(N - 3) + previousC(N - 1))/(n(N - 2));
        newC((N - 1):N) = (previousC((N - 3):(N - 2))/4 - previousC((N - 2):(N - 1)))./(n((N - 1):N));
        newC(N + 1) = previousC(N - 1)/(4*n(N + 1));
        newC = newC/4;
    end

end