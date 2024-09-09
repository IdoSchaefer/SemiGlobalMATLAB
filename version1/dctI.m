function dctv = dctI(v)
% The function performes a discrete cosine transform of the first kind.
% Suitable if it is desired to include the boundaries of the domain in
% the transform.
% The convention used for the factor of the sums, makes the transform its
% own inverse.
    dim = size(v);
    if dim(1) == 1
        N = dim(2);
        v = v.';
    else
        N = dim(1);
    end
    % The dctI is computed as a FFT of an extended vector, where v is
    % reflected at the right boundary (symmetric extension):
    dctv = 1/sqrt(2*(N - 1))*fft([v; v((N - 1):-1:2)]);
    dctv = dctv(1:N);
    if dim(1) == 1
        dctv = dctv.';
    end
    if isreal(v)
        dctv = real(dctv);
    elseif isreal(1i*v)
        dctv = 1i*imag(dctv);
    end
end