function mmiu = evmiu(psi, x, q)
    Nt = size(psi, 2);
    mmiu = zeros(1, Nt);
    for ti = 1:Nt
        mmiu(ti) = psi(:, ti)'*(x.*psi(:, ti));
    end
    if nargin>2
        mmiu = mmiu*q;
    end
end