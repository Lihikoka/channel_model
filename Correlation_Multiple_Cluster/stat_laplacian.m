function [composite_AoA_deg, composite_sigma_deg] = ... 
    stat_laplacian(number_clusters, Q, AoA_deg, sigma_deg, ... 
	delta_phi_deg) 
% [composite_AoA_deg, composite_sigma_deg] = ... 
%    stat_laplacian(number_clusters, Q, AoA_deg, sigma_deg, ... 
%    delta_phi_deg) 
% 
% Computes the composite AoA and sigma of the multi-cluster 
% Laplacian PAS 
% 
% 
% STANDARD DISCLAIMER 
% 
% CSys is furnishing this item "as is". CSys does not provide any 
% warranty of the item whatsoever, whether express, implied, or 
% statutory, including, but not limited to, any warranty of 
% merchantability or fitness for a particular purpose or any 
% warranty that the contents of the item will be error-free. 
% 
% In no respect shall CSys incur any liability for any damages, 
% including, but limited to, direct, indirect, special, or 
% consequential damages arising out of, resulting from, or any way 
% connected to the use of the item, whether or not based upon 
% warranty, contract, tort, or otherwise; whether or not injury was 
% sustained by persons or property or otherwise; and whether or not 
% loss was sustained from, or arose out of, the results of, the 
% item, or any services that may be provided by CSys. 
% 
% (c) Laurent Schumacher, AAU-TKN/IES/KOM/CPK/CSys - September 2002 
 
AoA_rad       = AoA_deg*pi/180; 
sigma_rad     = sigma_deg*pi/180; 
delta_phi_rad = delta_phi_deg*pi/180; 
 
composite_AoA_rad = sum(AoA_rad.*Q.*... 
                    (1-exp((-1).*sqrt(2).*(delta_phi_rad./sigma_rad)))); 
 
composite_sqr_sigma_rad = 0; 
for (k = 1:1:number_clusters) 
    z = sqrt(2)*(AoA_rad(k)-composite_AoA_rad)/sigma_rad(k); 
    composite_sqr_sigma_rad = composite_sqr_sigma_rad + ... 
                              .25*Q(k)*(sigma_rad(k)^2)*(z^2-2*z+2); 
    z = sqrt(2)*(AoA_rad(k)-composite_AoA_rad-delta_phi_rad(k))/sigma_rad(k); 
    composite_sqr_sigma_rad = composite_sqr_sigma_rad - ... 
                              .25*Q(k)*(sigma_rad(k)^2)*... 
                              exp((-1)*sqrt(2)*(delta_phi_rad(k)/sigma_rad(k)))*... 
                              (z^2-2*z+2); 
    z = (-1)*sqrt(2)*(AoA_rad(k)-composite_AoA_rad+delta_phi_rad(k))/sigma_rad(k); 
    composite_sqr_sigma_rad = composite_sqr_sigma_rad - ... 
                              .25*Q(k)*(sigma_rad(k)^2)*... 
                              exp((-1)*sqrt(2)*(delta_phi_rad(k)/sigma_rad(k)))*... 
                              (z^2-2*z+2); 
    z = (-1)*sqrt(2)*(AoA_rad(k)-composite_AoA_rad)/sigma_rad(k); 
    composite_sqr_sigma_rad = composite_sqr_sigma_rad + ... 
                              .25*Q(k)*(sigma_rad(k)^2)*(z^2-2*z+2); 
end; 
composite_AoA_deg = composite_AoA_rad*180/pi; 
composite_sigma_deg = sqrt(composite_sqr_sigma_rad)*180/pi;
