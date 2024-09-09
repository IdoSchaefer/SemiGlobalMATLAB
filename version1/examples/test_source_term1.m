% The program tests the efficiency of Hillel Tal-Ezer's
% semi-global propagator, for a forced harmonic oscillator with an arbitrary 
% inhomogeneous source term.
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
% An arbitrary inhomogeneous function:
ihfun = @(t) exp(-x.^2/2)*cos(t);
% The output time-grid:
dt = 0.1;
t=0:dt:T;
options = SGdefault_op;
data = SGdata(options);
% tic
% [U, mniter, matvecs, est_errors, history] = SemiGlobal1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, Nts, Nt_ts, Ncheb, tol, options, data);
% toc
display('Semi-global algorithm computation:')
tic
[U, mniter, matvecs, est_errors, history] = SemiGlobal1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, ihfun, [-195*1i, 0], fi0, t, Nts, Nt_ts, Ncheb, tol, options, data);
toc
% The mean number of iterations, for a time step (should be close to 1, for ideal
% efficiency):
mniter
% The number of matrix-vector multiplications, which constitutes the most of the
% numerical effort:
matvecs
% Checking if everything is OK - comparing to the ode45 result:
RKfun = @(t, u) -1i*Hpsi(K, V + x*cos(t), u) + ihfun(t);
options = odeset('RelTol', 1e-11, 'absTol', 1e-12);
display('ode45 computation:')
tic
[time, URK] = ode45(RKfun, [0 T/2 T], fi0, options);
toc
URK = URK(end, :).';
error = norm(U(:, end) - URK(:, end))/norm(URK(:, end))