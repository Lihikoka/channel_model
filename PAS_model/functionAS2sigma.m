function [sigmaA_L_in_degree, sigmaA_G_in_degree, sigma_L_in_degree, sigma_G_in_degree] = functionAS2sigma(sigma_in_degree)
    delta_phi_L = pi;
    delta_phi_G = pi;

    bound = 1000;
    
    phi_in_degree = -bound:0.01:bound;
    phi = phi_in_degree * pi / 180;

    sigma = sigma_in_degree * pi / 180;

    P_L = zeros(1, length(phi));
    P_G = zeros(1, length(phi));
    
    P_L_general = zeros(1, length(phi));
    P_G_general = zeros(1, length(phi));
    
    Q_L = 1/(1-exp(-sqrt(2)*delta_phi_L/sigma)); 
    Q_G = 1/erf(delta_phi_G/(sqrt(2)*sigma));

    sigmaA_L = 0;
    sigmaA_G = 0;
    sigma_L = 0;
    sigma_G = 0;
    cdf_G = 0;
    cdf_L = 0;


    for i=-bound:0.01:bound
        ii = i*pi/180;
        j = round((i+bound)*100+1);
        P_L(j) = Q_L/(sqrt(2)*sigma)*exp(-sqrt(2)*abs(ii)/sigma);
        P_G(j) = Q_G/(sqrt(2*pi)*sigma)*exp(-(ii)^2/(2*sigma^2));
        
        P_L_general(j) = 1/(sqrt(2)*sigma)*exp(-sqrt(2)*abs(ii)/sigma);
        P_G_general(j) = 1/(sqrt(2)*sigma)*exp(-sqrt(2)*abs(ii)/sigma);

        cdf_G = cdf_G + P_G_general(j)*(0.01*pi/180);
        cdf_L = cdf_L + P_L_general(j)*(0.01*pi/180);

        sigma_G = sigma_G + ii*ii*P_G_general(j)*(0.01*pi/180);
        sigma_L = sigma_L + ii*ii*P_L_general(j)*(0.01*pi/180);

        if i >= -delta_phi_L*180/pi && i <= delta_phi_L*180/pi
            sigmaA_L = sigmaA_L + ii*ii*P_L(j)*(0.01*pi/180);
        end
        
        if i >= -delta_phi_G*180/pi && i <= delta_phi_G*180/pi
            sigmaA_G = sigmaA_G + ii*ii*P_G(j)*(0.01*pi/180);
        end
        
    end

    sigma_L_in_degree = sqrt(sigma_L)*180/pi
    sigma_G_in_degree = sqrt(sigma_G)*180/pi

    sigma_in_degree

%     sqrt(sigma_L)
%     sqrt(sigma_G)
%     sigma
    
    sigmaA_L_in_degree = sqrt(sigmaA_L) * 180 / pi;
    sigmaA_G_in_degree = sqrt(sigmaA_G) * 180 / pi;
    
    if(sigma_in_degree == 70)
        figure(1)
        plot(phi_in_degree, P_L)
        hold on;
        plot(phi_in_degree, P_G)
        axis([-200, 200, 0, 5])
        grid on
    end
end