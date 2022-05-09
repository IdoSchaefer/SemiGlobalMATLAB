function u = RK4uf(fderiv, t0tf, u0, dt, varargin)
% The function solves a system of ode's by RK4.
% fderiv:  a function handle. The value of the function is the derivative of u with respect to t.
% t0tf:  a vector that contains initial and final time: [t0 tf]
% u0:  a vector of the initial values of the result vector u.
% dt: The time step.
% u:  The result u vector at the final time.
    t0 = t0tf(1);
    tf = t0tf(2);
    Nt = round(abs((tf - t0)/dt));
    u = u0;
    t = t0;
    for ti = 1:Nt
        ut = u;
 % The k's are estimations of du. The result is composed of a weighted
 % average of the 4 k's.
        k1 = feval(fderiv, t, ut, varargin{:})*dt;
        k2 = feval(fderiv, t + 0.5*dt, ut + 0.5*k1, varargin{:})*dt;
        k3 = feval(fderiv, t + 0.5*dt, ut + 0.5*k2, varargin{:})*dt;
        k4 = feval(fderiv, t + dt, ut + k3, varargin{:})*dt;
        u = u + (k1 + 2*k2 + 2*k3 + k4)/6;
        t = t + dt;
    end
end