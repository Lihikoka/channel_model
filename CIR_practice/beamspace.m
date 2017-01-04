clc
clear

Nt = 81; Nr = 1;  % Number of antenna at Tx and Rx
K = 20; % Number of user
L = 5; % Number of path
PL = [1, 0.3, 0.2, 0.15, 0.1]; % Path losses

phi = (rand(K, L)-0.5)*pi; % K users' azimuth in radian in [-pi/2, pi/2]
phi_degree = phi*180/pi;
theta = 0.5*sin(phi);
U = zeros(Nt, Nt);
I = zeros(Nt, 1);
a = zeros(Nt, L);
beta = zeros(L, 1);

% Beam index set
for l=1:Nt
    I(l) = (l-1) - (Nt-1)/2;
end

% Channel matrix
H_LOS = [];
H_NLOS = [];
H = [];

for k = 1:K
    for l = 1:L
        beta(l) = PL(l)*exp(-1i*2*pi*rand(1));
    end
    for n = 1:Nt
        for l = 1:L
            a(n, l) = exp(-1i*2*pi*theta(k,l)*I(n));
        end
    end
    h_k_LOS = beta(1)*a(:, 1);
    H_LOS = horzcat(H_LOS, h_k_LOS);
    h_k_NLOS = zeros(Nt, 1);
    for i = 2:L
        h_k_NLOS = h_k_NLOS + beta(i)*a(:, i);
    end
    H_NLOS = horzcat(H_NLOS, h_k_NLOS);
    H = H_LOS + H_NLOS;
end

% DFT matrix
for i = 1:Nt
    for j = 1:Nt
        U(i,j) = exp(-1i*2*pi*I(i)*I(j)/Nt)/sqrt(Nt);
    end
end

% Beamspace channel
H_b_LOS = ctranspose(U)*H_LOS; % LOS only
H_b_LOS_power = abs(ctranspose(H_b_LOS)).^2;

H_b = ctranspose(U)*H; % LOS and NLOS
H_b_power = abs(ctranspose(H_b)).^2;

H_b_NLOS = ctranspose(U)*H_NLOS; % NLOS only
H_b_NLOS_power = abs(ctranspose(H_b_NLOS)).^2;


figure(1)
contour(I*180/(Nt-1), 1:K, H_b_LOS_power);
title('Counter plot - LOS only');
xlabel('TX BEAM DIRECTION (degree)');
ylabel('MOBILE STATION INDEX');
axis([-90, 90 -inf, inf]);

figure(2)
contour(I*180/(Nt-1), 1:K, H_b_power);
title('Counter plot - LOS and NLOS');
xlabel('TX BEAM DIRECTION (degree)');
ylabel('MOBILE STATION INDEX');
axis([-90, 90 -inf, inf]);

figure(3)
contour(I*180/(Nt-1), 1:K, H_b_NLOS_power);
title('Counter plot - NLOS only');
xlabel('TX BEAM DIRECTION (degree)');
ylabel('MOBILE STATION INDEX');
axis([-90, 90 -inf, inf]);


for i=1:K
    figure(4)
    subplot(3, 1, 1);
    stem(I*180/(Nt-1), H_b_LOS_power(i, :));
    title('LOS only');
    xlabel('TX BEAM DIRECTION (degree)');
    ylabel('|H_b^H|^2');
    %axis([-90, 90 -inf, inf]);
    phi_degree(i,1)
    
    subplot(3, 1, 2);
    stem(I, H_b_power(i, :));
    title('LOS and NLOS');
    xlabel('TX BEAM INDEX');
    ylabel('|H_b^H|^2');
    
    subplot(3, 1, 3);
    stem(I, H_b_NLOS_power(i, :));
    title('NLOS only');
    xlabel('TX BEAM INDEX');
    ylabel('|H_b^H|^2');
    hold off
end