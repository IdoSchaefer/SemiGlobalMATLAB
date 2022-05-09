function Q = r2Taylor4(x, Dsize)
% The coefficients contain the 1/n! factor.
% The domain is of size 4 (capacity 1).
% x: A vector containing the sampling points
% Dsize: size of the domain.
    Dfactor = 4/Dsize;
    NC = length(x);
    Q = zeros(NC);
    % First row:
    Q(1, 1) = 1;
    % All the rest of the rows:
    for ri = 2:NC
        Q(ri, 1) = -Dfactor*x(ri - 1)*Q(ri - 1, 1);
        for Taylori = 2:ri-1
            Q(ri, Taylori) = (Q(ri - 1, Taylori - 1) - x(ri - 1)*Q(ri - 1, Taylori))*Dfactor;
        end
        Q(ri, ri) = Q(ri - 1, ri - 1)*Dfactor;
    end
end