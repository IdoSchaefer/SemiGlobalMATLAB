function allT = cheb_pols(x, N)
% The function computes the Chebyshev polynomials up to order N, evaluated
% at values x. The computation is performed by the recursive definition of
% the Chebyshev polynomials.
% x is a row vector.
   allT = zeros(N + 1, length(x));
   allT(1, :) = 1;
   allT(2, :) = x;
   for Ti = 3:(N + 1)
       allT(Ti, :) = 2*x.*allT(Ti - 1, :) - allT(Ti - 2, :);
   end
end