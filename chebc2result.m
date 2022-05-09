function result = chebc2result(Ccheb, xdomain, xresult)
% The function computes the Chebyshev polynomial approximation of a function
% from the corresponding Chebyshev coefficients, at a given set of points.
% Ccheb: The Chebyshev coefficients of the function (see chebc.m)
% xdomain: The approximation domain; a vector of the form [xmin, xmax]
% xresult: The set of points in which the function is to be evaluated
    % Transforming xresult to the Chebyshev domain [-1 1]:
    xrcheb = (2*xresult - xdomain(1) - xdomain(2))/(xdomain(2) - xdomain(1));
    m = length(Ccheb);
    % Pk represents the Chebyshev polynomial of the (k-1)'th degree.
    % P1 represents the Chebyshev polynomial of the (k-3)'th degree.
    % P2 represents the Chebyshev polynomial of the (k-2)'th degree.
    P1 = ones(size(xresult));
    P2 = xrcheb;
    result = Ccheb(1)*P1 + Ccheb(2)*P2;
    for k = 3:m
          Pk = 2*xrcheb.*P2 - P1;
          result = result + Ccheb(k).*Pk;
          P1 = P2;
          P2 = Pk;
    end
end