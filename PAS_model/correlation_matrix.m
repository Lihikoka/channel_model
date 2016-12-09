N = 4;
spacing = 0.5;
d = zeros(N);
for i=1:N
    for j=1:N
        d(i, j) = 0 + abs(i-j) * spacing;
    end 
end
phi_0 = 67.5;
AS = 35;
type = 'uniform';



R_MS = functionCorrelationCoefficient(d, phi_0, AS, type)