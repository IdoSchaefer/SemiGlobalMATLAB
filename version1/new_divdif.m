function [polcoef, diagonal_out] = new_divdif(allz, new_fz, diagonal_in)
% The function computes additional divided difference coefficients for a Newton
% interpolation, where several divided differences have already been
% computed. Applies when in addition to an existing set of sampling points,
% several new samplings are given. The new divided differences are computed
% from the new samplings and an existing diagonal in the divided difference
% table.
% Input:
% allz: A vector which contains the all the sampling points, including the
% old ones.
% new_fz: The function values at the new sampling points. Different
% sampling points are represented by different columns.
% diagonal_in: The last diagonal of the divided difference table, used for
% computation of the last old divided difference.
% Output:
% polcoef: The coefficients of the Newton basis polynomials for the Newton
% interpolation. The coefficients for the different Newton basis polynomials
% are represented by different columns.
% diagonal_out: The last diagonal, for continuation of the process to higher orders,
% if necessary.
    [dim, Nnewpoints] = size(new_fz);
    Noldpoints = size(diagonal_in, 2);
    Npoints = Noldpoints + Nnewpoints;
    polcoef = zeros(dim, Nnewpoints);
    diagonal_out = [diagonal_in, zeros(dim, Nnewpoints)];
    for coefi = (Noldpoints + 1):Npoints
        diagonal_out(:, coefi) = new_fz(:, coefi - Noldpoints);
        for dtermi = coefi-1:-1:1
            diagonal_out(:, dtermi) = (diagonal_out(:, dtermi + 1) - diagonal_out(:, dtermi))./(allz(coefi) - allz(dtermi));
        end
        polcoef(:, coefi - Noldpoints) = diagonal_out(:, 1);
    end
end