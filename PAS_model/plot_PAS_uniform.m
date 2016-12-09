clc
clear

AS_in_degree = [10 30];
AoA_in_degree = 0;

%%%%%%% Initialization %%%%%%%%
AS = zeros(length(AS_in_degree));
delta_phi = zeros(length(AS_in_degree));
Q = zeros(length(AS_in_degree));
delta_phi_in_degree = zeros(length(AS_in_degree));

for ii=1:length(AS_in_degree)
    AS(ii) = AS_in_degree(ii) * pi / 180;
    delta_phi(ii) = sqrt(3)*AS(ii);
    AoA = AoA_in_degree * pi / 180;
    Q(ii) = 1/(2*delta_phi(ii));
    delta_phi_in_degree(ii) = delta_phi(ii) * 180 / pi;

    phi_in_degree = -delta_phi_in_degree(ii)+AoA_in_degree : 0.1 : delta_phi_in_degree(ii)+AoA_in_degree;
    phi = -delta_phi(ii):(0.1*pi/180):delta_phi(ii);
    P = zeros(1, 200*2*10+1);
    
    sigmaA = 0;
    sigma = 0;
    j = 1;
    k = 1;
    for i=-200:0.1:200
        j = round((i+200)*10+1); % change phi_in_degree to index
        if(i < delta_phi_in_degree(ii) + AoA_in_degree && i > -delta_phi_in_degree(ii) + AoA_in_degree)
            P(j) = Q(ii);
            sigmaA = sigmaA + phi(k)^2*P(j)*(0.01*pi/180);
            sigma = sigma + phi(k)^2*P(j)*(0.01*pi/180);
            k = k + 1;
        else
            P(j) = 0;
        end
        
        j = j+1;
    end
    sigmaA
    sigmaA_in_degree = sqrt(sigmaA)*180/pi
   
    
    plot(-200:0.1:200, P);
    hold on
end