function C = recurrence_coefs(Cinit, Norders, f_rec, order_init)
% The function computes a set of coefficients which are defined by a
% recurrence relation, from a set of initial coefficients Cinit and a
% function which defines the recurrence rule.
% Input:
% Cinit: A matrix; contains a set of initial coefficients which are
% involved in the iterative process. The row index represents different
% orders of the recursive process, and the column index represents the
% different coeffcients of the same order in the recursive process.
% Note: In order to avoid repeated memory allocation, the number of columns
% should be the same as the final output matrix C.
% Norders: The number of orders in the recursive process to be computed
% f_rec: A function handle of the form @(previousC, prev_order);
% computes a new order of coefficients from the previous orders.
%   Input of f_rec:
%   previousC: A matrix of the previous coefficients of the form of Cinit.
%   prev_order: The order of previousC
%   Output of f_rec: a row vector of the new coefficients.
% Output:
% C: A matrix; the set of all coeffcients, including Cinit; the row and
% column index represent the same as in Cinit.
    if nargin<4
        order_init = 0;
    end
    [Nrec, Ncoefs] = size(Cinit);
    C = zeros(Norders + Nrec, Ncoefs);
    C(1:Nrec, :) = Cinit;
    for Ci = (Nrec + 1):(Nrec + Norders)
        C(Ci, :) = f_rec(C((Ci - Nrec):(Ci - 1), :), order_init + Ci - 1);
    end
end