[Uex, mniterex, matvecsex] = SemiGlobalArnoldi_xp(K240, Vabs240, @(u,x,t) -xabs240*0.1*sech((t-500)/(170)).^2.*cos(0.06*(t-500)), [], fi0240, x240, [0, 1000], 3e4, 9, 13, eps);