function [V, H] = createKrop(Op, v0, Nv)
% The function creates an orthonormal Krylov basis of dimension Nv+1, from 
% a vector v0, and a function handle Op, which represents the operation of
% an operator on a vector. It uses the Arnoldi algorithm.
% The columns of V are the vectors of the Krylov space.
% H is the extended Hessenberg matrix.
    lenv = length(v0);
    V = zeros(lenv, Nv+1);
    H = zeros(Nv+1, Nv);
    V(:, 1) = v0/norm(v0);
    for vj = 1:Nv
        V(:, vj+1) = Op(V(:, vj));
        for vi = 1:vj
            H(vi, vj) = V(:, vi)'*V(:, vj+1);
            V(:, vj+1) = V(:, vj+1) - H(vi, vj)*V(:, vi);
        end
        H(vj+1, vj) = norm(V(:, vj+1));
        V(:, vj+1) = V(:, vj+1)/H(vj+1, vj);
    end
end