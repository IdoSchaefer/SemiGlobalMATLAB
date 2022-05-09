function M = vchebMop(operator, u0, leftb, rightb, Ncheb)
% The function computes the vectors which can be used for a Chebyshev
% expansion of any function of an operator that operates on the vector u0.
% The expansion will be: f(operator)u0 = sum(a_i*M(:, i+1)), i = 0, 1, ..., Ncheb - 1 
% The coefficients a_i are f dependent, and will be computed by the
% program chebc.
% The program is useful when we want to compute a number of functions
% of the same operator, which operate on the same vector.
% operator: a function handle of the form: @(v). The function will return 
% the operation of the operator on the vector v.
% The eigenvalue domain of the operator is: [leftb rightb].
    % Defining the operator in the domain of the Chebyshev polynomials:
    % [-1 1]
    chebop = @(v) (2*operator(v) - (leftb + rightb)*v)/(rightb - leftb);
    dim = length(u0);
    M = zeros(dim, Ncheb);
    M(:, 1) = u0;
    M(:, 2) = chebop(u0);
    for k = 3:Ncheb
        M(:, k) = 2*chebop(M(:, k-1)) - M(:, k-2);
    end
end