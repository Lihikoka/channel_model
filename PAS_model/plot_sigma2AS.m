clc
clear

sigma_in_degree = 0:1:120;

sigmaA_L_in_degree = zeros(1, length(sigma_in_degree));
sigmaA_G_in_degree = zeros(1, length(sigma_in_degree));
sigma_L_in_degree = zeros(1, length(sigma_in_degree));
sigma_G_in_degree = zeros(1, length(sigma_in_degree));

for i=1:length(sigma_in_degree)
    [sigmaA_L_in_degree(i), sigmaA_G_in_degree(i), sigma_L_in_degree(i), sigma_G_in_degree(i)] = functionAS2sigma(sigma_in_degree(i));
end
figure(2)
plot(sigmaA_G_in_degree, sigma_in_degree);
hold on
plot(sigmaA_L_in_degree, sigma_in_degree);
axis([0, 80, 0, 120])
grid on
legend('Gaussian PAS', 'Laplacian PAS');

