function [gs, niter] = gsNLHdiag(H0, Vnl, x, tol)
% The function finds the ground state of a non-linear Hamiltonian, by an
% iterative process.
% H0 is the linear part of the Hamiltonian, represented by a matrix.
% Vnl is the nonlinear purterbation. It's a function handle of the form:
% @(u, x, t).
% x is the space grid.
% tol is the desired tolerance.
    szx = size(x);
    if szx(1) == 1
        x = x.';
    end
    [P, D] = eig(H0);
    eigval = diag(D);
    [minE, nminE] = min(eigval);
    gs = P(:, nminE);
    H = H0 + diag(Vnl(gs, x, 0));
    niter = 1;
    maxNiter = 100;
    while abs(exval(H^2, gs) - exval(H, gs)^2)>tol && niter<=maxNiter
        [P, D] = eig(H);
        eigval = diag(D);
        [minE, nminE] = min(eigval);
        gs = P(:, nminE);
        H = H0 + diag(Vnl(gs, x, 0));
        niter = niter + 1;
    end
    if niter>maxNiter
        fprintf('The program has failed to achieve the desired tolerance.\n');
    end
    niter
end