funciton [R R_BS R_MS] = geometry2correlation()
% geometry2correlation
%
% User interface to the definition of the correlation matrices
% of the MIMO radio channel model. For both Node B and UE,
% the user is prompted for
%
% * the number of elements
% * the spacing between elements (Uniform Linear Array is assumed)
% * the number of impinging clusters of waves
% * their Power Azimuth Spectrum (PAS) type (Uniform, Gaussian or Laplacian)
% * the mean angle of incidence (in degrees)
% * the Azimuth Spread (AS)
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
% (c) Laurent Schumacher, AAU/TKN/I8/KOM/CPK/CSys - September 2001

clear;
clf;
display = 1;

%
% Interactive dialog
%

% Downlink or uplink

string = sprintf('\n');
disp(string);
Downlink = -1;
while ((Downlink ~= 0) && (Downlink ~=1))
    Downlink = input('Downlink (Yes = 1, no = 0)? ');
end;

% Field or envelope/power

string = sprintf('\n');
disp(string);
Field = -1;
while ((Field ~= 0) && (Field ~=1))
    Field = input('Type of correlation coeffi1cient (Field = 0, Envelope = 1)? ');
end;

% Configuration of Node B

string = sprintf('\n*** Node B ***\n');
[M, spacing_Node_B, d_norm_Node_B, cluster_number_Node_B, ...
 amplitude_cluster_Node_B, PAS_type_Node_B, phi_deg_Node_B, ...
 AS_deg_Node_B, delta_phi_deg_Node_B] = dialog(string);

% Configuration of UE

string = sprintf('\n*** UE ***\n');
[N, spacing_UE, d_norm_UE, cluster_number_UE, ...
 amplitude_cluster_UE, PAS_type_UE, phi_deg_UE, AS_deg_UE, ...
 delta_phi_deg_UE] = dialog(string);

%
% Computation of the correlation matrices
%

% Node B

[R_Node_B, Q_Node_B, sigma_deg_Node_B] = correlation(M, spacing_Node_B, ...
						 d_norm_Node_B, ...
						 cluster_number_Node_B, ...
						 amplitude_cluster_Node_B, ...
						 PAS_type_Node_B, ...
						 phi_deg_Node_B, ...
						 AS_deg_Node_B, ...
                         delta_phi_deg_Node_B, Field);

if (display && (~(M==1)))
  figure(1);
  clf;
  switch(PAS_type_Node_B)
   case -1
    plot_uniform(cluster_number_Node_B, Q_Node_B, phi_deg_Node_B, ...
		 AS_deg_Node_B, 1, 'b-');
   case 2
    plot_gaussian(cluster_number_Node_B, Q_Node_B, phi_deg_Node_B, ...
		  sigma_deg_Node_B, delta_phi_deg_Node_B, 1, 'b-');
   case 3
    plot_laplacian(cluster_number_Node_B, Q_Node_B, phi_deg_Node_B, ...
		   sigma_deg_Node_B, delta_phi_deg_Node_B, 1, 'b-');
%    otherwise
%     break;
  end;
end;

% UE

[R_UE, Q_UE, sigma_deg_UE] = correlation(N, spacing_UE, d_norm_UE, ...
					 cluster_number_UE, ...
					 amplitude_cluster_UE, ...
					 PAS_type_UE, phi_deg_UE, ...
                     AS_deg_UE, delta_phi_deg_UE, Field);
if (display && (~(N==1)))
  figure(1);
  switch(PAS_type_UE)
   case 1
    plot_uniform(cluster_number_UE, Q_UE, phi_deg_UE, AS_deg_UE, 2, 'b-');
   case 2
    plot_gaussian(cluster_number_UE, Q_UE, phi_deg_UE, sigma_deg_UE, ...
        delta_phi_deg_UE, 2, 'b-');
   case 3
    plot_laplacian(cluster_number_UE, Q_UE, phi_deg_UE, sigma_deg_UE, ...
		delta_phi_deg_UE, 2, 'b-');
%    default 
%     break
  end;
end;

%
% Information exploitation
%

if (Downlink)
    nTx = M;
    nRx = N;
    R   = kron(R_Node_B, R_UE);

else
    nTx = N;
    nRX = M;
    R   = kron(R_UE, R_Node_B);
end;