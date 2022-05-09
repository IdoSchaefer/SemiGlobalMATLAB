function result = exval(A, psi)
% The function returns the expectation value of the operator A, for an unnormalized wave
% function, psi.
% A is a matrix, and psi is a vector.
    result = (psi'*A*psi)/(psi'*psi);
end

