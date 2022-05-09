% The program tests the efficiency of Hillel Tal-Ezer's
% semi-global propagator, for a forced harmonic oscillator.
% You may play with the following parameters:
T = 10; Nts = 200; Nt_ts = 9; Ncheb = 9; tol=1e-5;
% Constructing the grid:
L = 16*sqrt(pi);
Nx = 128;
dx = L/Nx;
x = (-L/2:dx:(L/2 - dx)).';
% Constructing the kinetic energy matrix diagonal in the p domain:
p = (0:(2*pi/L):(2*pi*(1/dx - 1/L))).';
p((Nx/2 + 1):Nx) = p((Nx/2 + 1):Nx) - 2*pi/dx;
K = p.^2/2;
% The potential energy matrix diagonal in the x domain:
V = x.^2/2;
% The harmonic oscillator ground state:
fi0 = pi^(-1/4)*exp(-x.^2/2)*sqrt(dx);
% The output time-grid:
dt = 0.1;
t=0:dt:T;
tic
[U, mniter, matvecs] = SemiGlobal(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - ones(1, Nt_ts)*cos(t2)).*u1, 0, [], [-188*1i, 1i], fi0, t, Nts, Nt_ts, Ncheb, tol);
toc
% The mean number of iterations, for a time step (should be close to 1, for ideal
% efficiency):
mniter
% The number of matrix-vector multiplications, which constitutes the most of the
% numerical effort:
matvecs
% Computation of the maximal error - the deviation from the analytical
% result of the expectation value.
% Computing the expectation value of x in all the time points:
mx = evmiu(U, x);
error = mx - (-0.5*sin(t).*t);
maxer = max(abs(error))