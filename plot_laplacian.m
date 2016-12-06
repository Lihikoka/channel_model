function plot_laplacian(number_clusters, power, AOA_deg, sigma_deg, ... 
		      delta_phi_deg, plot_number, colour) 
% plot_laplacian(number_clusters, power, AOA_deg, sigma_deg, 
%                delta_phi_deg, plot_number, colour) 
% 
% Plots the multi-cluster Laplacian PAS 
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
 
phi_deg   = -180:0.5:180; 
maximum   = 0; 
 
phi_rad   = phi_deg*pi/180; 
sigma_rad = sigma_deg*pi/180; 
 
for (k = 1:1:number_clusters) 
  mask = .5.*(sign(delta_phi_deg(k)-phi_deg+1e-6)+... 
              sign(delta_phi_deg(k)+phi_deg+1e-6)); 
  PAS = (power(k)/(sqrt(2)*sigma_rad(k))).*exp(((-1).*sqrt(2).*abs(phi_rad))./sigma_rad(k)); 
  if (AOA_deg(k) > 0) 
    limit = find(phi_deg==(180-AOA_deg(k))); 
    mask  = [mask(limit+1:size(mask,2)),mask(1:limit)]; 
    PAS = [PAS(limit+1:size(PAS,2)),PAS(1:limit)]; 
  else 
    limit = find(phi_deg==(-180-AOA_deg(k))); 
    mask  = [mask(limit:size(mask,2)),mask(1:limit-1)]; 
    PAS = [PAS(limit:size(PAS,2)),PAS(1:limit-1)]; 
  end; 
  maximum = max(maximum,max(PAS)); 
  subplot(2,1,plot_number),plot(phi_deg, mask.*PAS, colour); 
  hold on; 
end; 
axis([-180 180 0 1.1*maximum]); 
grid; 
ylabel('Power (Linear)'); 
if (plot_number == 1) 
  title('Power Azimuth Spectrum'); 
else 
  xlabel('\phi [degrees]'); 
end;
