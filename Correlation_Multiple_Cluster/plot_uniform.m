function plot_uniform(number_clusters, amplitude, AOA, AS, ... 
		      plot_number, colour) 
% plot_uniform(number_clusters, amplitude, AOA, AS, ... 
%              plot_number, colour) 
% 
% Plots the multi-cluster uniform PAS 
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
 
phi_deg = -180:1:180; 
for (k = 1:1:number_clusters) 
  cluster_min = AOA(k) - sqrt(3)*AS(k); 
  cluster_MAX = AOA(k) + sqrt(3)*AS(k); 
  if (cluster_min < -180) 
    span = [-180:1:AOA(k) + sqrt(3)*AS(k),360+AOA(k) - sqrt(3)* ... 
	    AS(k):1:180]; 
    cluster_min = cluster_min + 360; 
  elseif (cluster_MAX > 180) 
    span = [-180:1:AOA(k) + sqrt(3)*AS(k)-360,AOA(k) - sqrt(3)* ... 
	    AS(k):1:180]; 
    cluster_MAX = cluster_MAX - 360; 
    else 
      span = AOA(k) - sqrt(3)*AS(k):1:AOA(k) + sqrt(3)*AS(k); 
  end; 
  subplot(2,1,plot_number),plot(span, amplitude(k), colour); 
  hold on; 
  subplot(2,1,plot_number),lh = line([cluster_min cluster_min], ... 
      [0 amplitude(k)]); 
  set(lh, 'Color', colour(1), 'LineStyle', colour(2:size(colour,2))); 
  subplot(2,1,plot_number),lh = line([AOA(k) AOA(k)], ... 
      [0 amplitude(k)]); 
  set(lh, 'Color', colour(1), 'LineStyle', colour(2:size(colour,2))); 
  subplot(2,1,plot_number),lh = line([cluster_MAX cluster_MAX], ... 
      [0 amplitude(k)]); 
  set(lh, 'Color', colour(1), 'LineStyle', colour(2:size(colour,2))); 
end; 
axis([-180 180 0 1.1*max(amplitude)]); 
grid; 
ylabel('Power (Linear)'); 
if (plot_number == 1) 
  title('Power Azimuth Spectrum'); 
else 
  xlabel('\phi [degrees]'); 
end;
