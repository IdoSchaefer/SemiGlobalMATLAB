% The program tests the efficiency of Hillel Tal-Ezer's
% semi-global propagator, for a forced harmonic oscillator in a BEC trap.
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
% We have to construct the H0 matrix, to find the ground state (it's unnecessary
% for the propagation process).
% The potential energy matrix:
Vmat = diag(V);
% The kinetic energy matrix in the x domain:
Kmat = Nx*ifft(ifft(diag(K))')';
% The Hamiltonian:
H = Kmat + Vmat;
% The ground state, found by an iterative process:
gs = gsNLHdiag(H, @(u,x,t) conj(u).*u, x, 2e-12);
% The output time-grid:
dt = 0.1;
t=0:dt:T;
options = SGdefault_op;
data = SGdata(options);
tic
[U, mniter, matvecs, est_errors, history] = SemiGlobal1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t) + conj(u).*u, v), @(u1, t1, u2, t2) -1i*(x*(cos(t1) - cos(t2)) + conj(u1).*u1 - (conj(u2).*u2)).*u1, 0, [], [-195*1i, 0], gs, t, Nts, Nt_ts, Ncheb, tol, options, data);
toc
% The mean number of iterations, for a time step (should be close to 1, for ideal
% efficiency):
mniter
% The number of matrix-vector multiplications, which constitutes the most of the
% numerical effort:
matvecs
fprintf('\nSmall tolerance Runge-Kutta computation:\n')
optionsRK = odeset('RelTol', 1e-13, 'absTol', 1e-13);
tic
[time, URK] = ode45(@(t, u) -1i*Hpsi(K, V + x*cos(t) + conj(u).*u, u), [0 T/2 T], gs, optionsRK);
toc
URK = URK(end,:).';
% Error from the RK solution:
error = norm(U(:, end) - URK(:, end))/norm(URK(:, end))