clc
clear

delta_phi = pi;
delta_phi_deg = delta_phi * 180 / pi;
bound = 200;
sigmaA_deg = [10 30];
phi_0 = 0;  % AoA


for iter=1:length(sigmaA_deg)
    string = ['load AS2sigma_laplacian_', num2str(delta_phi_deg), '.mat'];
    eval(string);

    phi_deg = -bound:0.01:bound;
    phi = phi_deg * pi / 180;

    index = find(round(sigmaA_L_deg)==sigmaA_deg(iter));
    sigma = sigma_deg(index) * pi / 180;
    P_L = zeros(1, length(phi));
    %P_L_general = zeros(1, length(phi));
    Q_L = 1/(1-exp(-sqrt(2)*delta_phi/sigma)); 

    sigmaA_L = 0;
    sigma_L = 0;
    cdf_L = 0;


    for i=-bound:0.01:bound
        ii = i*pi/180;
        j = round((i+bound)*100+1);
        P_L(j) = Q_L/(sqrt(2)*sigma)*exp(-sqrt(2)*abs(ii)/sigma);
        %P_L_general(j) = 1/(sqrt(2)*sigma)*exp(-sqrt(2)*abs(ii)/sigma);
        %cdf_L = cdf_L + P_L_general(j)*(0.01*pi/180);

        if i >= -delta_phi*180/pi && i <= delta_phi*180/pi
            sigmaA_L = sigmaA_L + ii*ii*P_L(j)*(0.01*pi/180);
        end
    end

    Rxx = zeros(1, length(0.01:0.01:6));
    for d=0.01:0.01:6
        Rxx(1) = besselj(0, d);
        i = round(d * 100);
        for m=1:100
            Rxx(i) = Rxx(i) + 4 * Q_L * besselj(2*m, d) * cos((2*m)*phi_0) * ...
                        (sqrt(2)/sigma + exp(-sqrt(2)*delta_phi/sigma) * (2*m*sin(2*m*delta_phi)-sqrt(2)*cos(2*m*delta_phi)/sigma)) / ...
                        (sqrt(2)*sigma*((sqrt(2)/sigma)^2+(2*m)^2));
        end
    end
    
    figure(1)
    plot(phi_deg, P_L)
    hold on;
    
    figure(2)
    plot(0.01:0.01:6, abs(Rxx));
    hold on;
end
figure(1)
axis([-200, 200, 0, 5])
grid on
legend('AS=10(degree)', 'AS=30(degree)');

figure(2)
grid on
legend('AS=10(degree)', 'AS=30(degree)');

