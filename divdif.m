function [polcoef, diagonal] = divdif(z, fz)
% The function computes the divided difference coefficients for a Newton interpolation.
% The routine is based on a divided difference table, where each
% coefficient is given by the last term of a new diagonal in the table.
% The program applies also for an interpolated function which returns a vector.
% Input:
% z: A vector which contains the sampling points
% fz: The function values at z. Different sampling points are represented
% by different columns.
% Output:
% polcoef: The coefficients of the Newton basis polynomials for the Newton
% interpolation. The coefficints for the different Newton basis polynomials
% are represented by different columns.
% diagonal: The last diagonal, for continuing the process to higher orders,
% if necessary.

    [dim, Npoints] = size(fz);
    polcoef = zeros(dim, Npoints);
    diagonal = zeros(dim, Npoints);
    polcoef(:, 1) = fz(:, 1);
    diagonal(:, 1) = fz(:, 1);
    % coefi indexes the Newton interpolation coefficients.
    % dtermi indexes the terms of the new diagonal in the divided difference
    % table. dtermi coincides with the index of the sampling point with the lowest
    % index of the divided difference. For example, the index of f[z(2), z(3), z(4)]
    % is dtermi = 2.
    for coefi = 2:Npoints
        diagonal(:, coefi) = fz(:, coefi);
        for dtermi = coefi-1:-1:1
            % diagonal(:, dtermi) belongs to the previous diagonal.
            % diagonal(:, dtermi + 1) belongs to the new diagonal.
            diagonal(:, dtermi) = (diagonal(:, dtermi + 1) - diagonal(:, dtermi))./(z(coefi) - z(dtermi));
        end
        polcoef(:, coefi) = diagonal(:, 1);
    end
end