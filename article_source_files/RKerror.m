T = 1e3;
Nsamp = 22;
minNt = 1.1e4;
allNtRK = zeros(Nsamp, 1);
allmvRK = zeros(Nsamp, 1);
allerRK = zeros(Nsamp, 1);
load coulomb_optV240
load Uex_article
for degi = 1:Nsamp
    deg = log10(minNt) + (degi-1)*0.1;
    Nt = round(10^deg);
    allNtRK(degi) = Nt;
    dt = T/Nt;
    uRK = RK4uf(@SEnl, [0 T], fi0240, dt, K240, @(u, x, t) Vabs240 -xabs240*0.1*sech((t-500)/(170)).^2.*cos(0.06*(t-500)), x240);
    matvecs = Nt*4;
    allmvRK(degi) = matvecs;
    if ~isfinite(uRK)
        display('Error.')
    end
    error = norm(uRK - Uex(:, end))/norm(Uex(:, end));
    allerRK(degi) = error;
end
figure
plot(log10(allmvRK), log10(allerRK), '-o')
xlabel('log(matvecs)')
ylabel('log(error)')
RKlincoef=polyfit(log10(allmvRK(1:21)), log10(allerRK(1:21)), 1);