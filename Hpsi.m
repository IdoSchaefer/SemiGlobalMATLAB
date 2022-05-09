function result = Hpsi(K, V, psi)
% The result is the operation of the Hamiltonian on the wave function
% psi.
% V is the potential energy. It's a vector in the x domain.
% K is the kinetic energy vector in the p domain. It has to be computed in the
% following way:
%     p = (0:(2*pi/xdlength):(2*pi*(1/dx - 1/xdlength))).';
%     p((Nx/2 + 1):Nx) = p((Nx/2 + 1):Nx) - 2*pi/dx;
%     K = p.^2/2;
%(Look in RKSE2)
    result = ifft(K.*fft(psi)) + V.*psi;
end