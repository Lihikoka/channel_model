function R = functionCorrelationCoefficient(d, phi_0, AS, type)

% d: an matrix that contains the distance between any two antennas at BS or MS
% phi_0: AoA in degree
% AS: Azumith Spread or sigmaA in degree
% type: which PAS model to use(uniform, truncated_gaussian or truncated_laplacian)


D = 2*pi*d; % normalized antenna distance
phi_0_radian = phi_0 * pi / 180;
sigmaA = AS * pi / 180;
Rxx = zeros(size(D));
Rxy = zeros(size(D));

switch lower(type)
    case 'uniform'
        delta_phi = sqrt(3)*sigmaA;
        Q = 1/(2*delta_phi);
        for m=1:100
            for i=1:size(D, 1)
                for j=1:size(D, 2)
                    Rxx(i, j) = Rxx(i, j) + besselj(0, D(i, j)) + 4 * Q * besselj(2*m, D(i, j)) * cos(2*m*phi_0_radian) * cos(2*m*delta_phi) / (2*m);
                    Rxy(i, j) = Rxy(i, j) + 4 * Q * besselj(2*m+1, D(i, j)) * sin((2*m+1)*phi_0_radian) * sin((2*m+1)*delta_phi) / (2*m+1);
                end
            end
        end
        R = Rxx + 1i * Rxy;
        
    case 'truncated_gaussian'
        delta_phi = sqrt(3)*sigmaA;
        Q = 1/(2*delta_phi);
        for m=1:100
            for i=1:size(D, 1)
                for j=1:size(D, 2)
                    Rxx(i, j) = Rxx(i, j) + besselj(0, D(i, j)) + 4 * Q * besselj(2*m, D(i, j)) * cos(2*m*phi_0) * cos(2*m*delta_phi) / (2*m);
                    Rxy(i, j) = Rxy(i, j) + 4 * Q * besselj(2*m+1, D(i, j)) * sin((2*m+1)*phi_0) * sin((2*m+1)*delta_phi) / (2*m+1);
                end
            end
        end
        R = Rxx + 1i * Rxy;
        
    case 'truncated_laplacian'
        delta_phi = pi;
        Q = 1 / erf(delta_phi / (sqrt(2)*sigma));
        for m=1:100
            for i=1:size(D, 1)
                for j=1:size(D, 2)
                    Rxx(i, j) = Rxx(i, j) + besselj(0, D(i, j)) +4 * Q * besselj(2*m, D(i, j)) * cos(2*m*phi_0) * cos(2*m*delta_phi) / (2*m);
                    Rxy(i, j) = Rxy(i, j) + 4 * Q * besselj(2*m+1, D(i, j)) * sin((2*m+1)*phi_0) * sin((2*m+1)*delta_phi) / (2*m+1);
                end
            end
        end
        R = Rxx + 1i * Rxy;
        
end
end