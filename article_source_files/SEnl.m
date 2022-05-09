function Du = SEnl(t, u, K, Vt, x)
    Du = -1i*(ifft(K.*fft(u)) + Vt(u, x, t).*u);
end