clc
clear

n = 81;
K = 1;
a = zeros(n, 1);
I = zeros(n, 1);
for l=1:n
    I(l) = (l-1) - (n-1)/2;
end

U = zeros(n, n);
for i = 1:n
    for j = 1:n
        U(i,j) = exp(-1i*2*pi*I(i)*I(j)/n)/sqrt(n);
    end
end

phi_degree = -90:90;
phi = phi_degree*pi/180;
theta = 0.5*sin(phi);
for ph = 1:length(phi_degree)
    for i = 1:n
        a(i) = exp(-1i*2*pi*theta(ph)*I(i));
    end
    H = a;

    H_b = ctranspose(U)*H;
    H_b_amp = abs(H_b);
    H_b_phase = log2(H_b);
    H_b_power = abs(ctranspose(H_b)).^2;
    figure(1)
    stem(asin(2*I/(n-1))*180/pi, H_b_power);
    phi_degree(ph)
   
    % figure(2)
    % stem(I, H_b_power);
end