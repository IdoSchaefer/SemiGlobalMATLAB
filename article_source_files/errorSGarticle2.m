function [allNt, allmv, aller, max_ers] = errorSGarticle2(T, Nt_ts, Nkr, minNt, Nsamp, Niter)
% The function computes the data of the error decay with reduction of the
% time-step.
if nargin<6
    Niter = 1;
end
allNt = zeros(Nsamp, 1);
allmv = zeros(Nsamp, 1);
aller = zeros(Nsamp, 1);
max_ers.texp = zeros(Nsamp, 1);
max_ers.FU = zeros(Nsamp, 1);
max_ers.conv = zeros(Nsamp, 1);
load coulomb_optV240
load Uex_article
for degi = 1:Nsamp
    deg = log10(minNt) + (degi-1)*0.1;
    Nt = round(10^deg);
    allNt(degi) = Nt;
    [u, ~, matvecs, all_est_er] = SemiGlobalArnoldi_xp(K240, Vabs240, @(u,x,t) -xabs240*0.1*sech((t-500)/(170)).^2.*cos(0.06*(t-500)), [], fi0240, x240, [0 T], Nt, Nt_ts, Nkr, eps, Niter, 20, false);
    allmv(degi) = matvecs;
    max_ers.texp(degi) = all_est_er.texp;
    max_ers.fU(degi) = all_est_er.fU;
    max_ers.conv(degi) = all_est_er.conv;
    if ~isfinite(u(:,2))
        display('Error.')
    end
    error = norm(u(:, 2) - Uex(:, end))/norm(Uex(:, end));
    aller(degi) = error;
end
figure
plot(log10(allmv), log10(aller), '-o')
xlabel('log(matvecs)')
ylabel('log(error)')
end