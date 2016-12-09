clc
clear

AS_in_degree = [10];
AoA_in_degree = 0;

AS = zeros(length(AS_in_degree));
delta_phi = zeros(length(AS_in_degree));
Q = zeros(length(AS_in_degree));
delta_phi_in_degree = zeros(length(AS_in_degree));

for ii=1:length(AS_in_degree)
    AS(ii) = AS_in_degree(ii) * pi / 180;
    delta_phi(ii) = pi;
    AoA = AoA_in_degree * pi / 180;
    Q(ii) = 1/erf(delta_phi(ii)/sqrt(2));
    delta_phi_in_degree(ii) = delta_phi(ii) * 180 / pi;

    phi_in_degree = -delta_phi_in_degree(ii):0.1:delta_phi_in_degree(ii);
    phi = -delta_phi(ii):(0.1*pi/180):delta_phi(ii);
    P = zeros(1, length(phi));
    
    sigmaA = 0;
    
    for i=1:length(phi)
        P(i) = Q(ii)/(sqrt(2*pi)*sigma)*exp(-(phi(i)-AoA)^2/(2*sigma^2));
        sigmaA = sigmaA + phi(i)^2*P(i)*(0.1*pi/180);
    end
    sigmaA
    sigmaA_in_degree = sqrt(sigmaA)*180/pi
   
    
    plot(phi_in_degree, P);
    hold on
end
%%%%%%%%%%%%%%%%%%%%%%
