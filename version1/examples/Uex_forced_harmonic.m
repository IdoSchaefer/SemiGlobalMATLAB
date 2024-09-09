% Analytic solution for the drive of test_harmonic.m
mx_ex = (-0.5*sin(t).*t);
mp_ex = -0.5*(sin(t) + t.*cos(t));
angle_analytic = (t/2 - (sin(2*t)/4 - t.*cos(2*t)/2)/8);
Uex_phase = pi^(-1/4)*exp(-1i*angle_analytic).*exp(1i*(mp_ex.*(x - mx_ex/2)) - (x - mx_ex).^2/2)*sqrt(dx);
Uex_fHO = Uex_phase(:,end);