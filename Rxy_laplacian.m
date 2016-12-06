function result = Rxy_laplacian(d_norm,phi_0_deg,sigma_deg,delta_phi_deg) 
% result = Rxy_laplacian(d_norm,phi_0_deg,sigma_deg,delta_phi_deg) 
% 
% Computes the correlation of the real and the imaginary part of 
% the signals in the case of a laplacian Power Azimuth Spectrum 
% (PAS) at spacings given by d_norm. The PAS is characterised by 
% the Angle Of Arrival (AOA) phi_0_deg and by its standard 
% deviation sigma_deg. 
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
% (c) Laurent Schumacher, AAU-TKN/IES/KOM/CPK/CSys - July 2001 
 
D             = 2*pi*d_norm; 
% Conversion degree -> rad 
phi_0_rad     = phi_0_deg*(pi/180); 
sigma_rad     = sigma_deg*(pi/180); 
delta_phi_rad = delta_phi_deg*(pi/180); 
 
m      = 0; 
result = (4.*besselj(1,D).*sin(phi_0_rad) ... 
	  .*((sqrt(2)./sigma_rad)+(exp((-1)*sqrt(2)*delta_phi_rad/sigma_rad) ... 
	        *(sin(delta_phi_rad)-sqrt(2)*cos(delta_phi_rad)/ ... 
		  sigma_rad)))) ... 
	 ./(sqrt(2)*sigma_rad*((sqrt(2)/sigma_rad)^2+1)); 
tmp_xy_laplacian = ones(size(result)); 
while (m < 100) 
    m                = m +1; 
    tmp_xy_laplacian =  (4.*besselj(2*m+1,D).*sin((2*m+1)*phi_0_rad) ... 
			 .*((sqrt(2)./sigma_rad)+(exp((-1)*sqrt(2)* ... 
						     delta_phi_rad/ ... 
						     sigma_rad) ... 
                               *((2*m+1)*sin((2*m+1)*delta_phi_rad)-sqrt(2)*cos((2*m+1)*delta_phi_rad)/sigma_rad)))) ... 
			./(sqrt(2)*sigma_rad*((sqrt(2)/sigma_rad)^2+(2*m+1)^2)); 
    result = result + tmp_xy_laplacian; 
end; 