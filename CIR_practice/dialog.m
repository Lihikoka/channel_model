function [M, spacing, d_norm, cluster_number, amplitude_cluster, ...
          PAS_type, phi_deg, AS_deg, delta_phi_deg] = dialog(string)
% [M, spacing, d_norm, cluster_number, amplitude_cluster,
% PAS_type, phi_deg, AS_deg, delta_phi_deg] = dialog(string)
%
% User interface to the definition of the MIMO radio channel
% model. The user is prompted for
%
% * the number of elements
% * the spacing between elements (Uniform Linear Array is assumed)
% * the number of impinging clusters of waves
% * their relative amplitude
% * their half domain definition (in degrees)
% * their Power Azimuth Spectrum (PAS) type (Uniform, Gaussian or
%   Laplacian)
% * their half domain definition (in degrees, usually 180 since the
%   PAS is defined over [ AOA-180 degrees, AOA + 180 degrees ]
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
% (c) Laurent Schumacher, AAU/TKN/I8/KOM/CPK/CSys - July 2001

disp(string);
M = 0;
while (M <1)
    M = input('Number of elements? ');
end;
if (M>1)
  spacing = 0;
  while (spacing <= 0)
    spacing = input('Normalised spacing between elements (in wavelengths)? ');
  end;
  d_norm = linspace(0,(M-1)*spacing,M);
  cluster_number = 0;
  while (cluster_number < 1)
    cluster_number = input(['Number of impinging clusters of waves? ']);
  end;
  PAS_type = 0;
  while ((PAS_type < 1) | (PAS_type > 3))
    PAS_type = input(['PAS type (1 = Uniform, 2 = Gaussian, 3 = Laplacian)? ']);
  end;
  if (cluster_number == 1)
    if ~(PAS_type == 1)
        delta_phi_deg = -1;
        while (delta_phi_deg <= 0)
            delta_phi_deg = input('Half Domain Definition (in degrees)? ');
        end;
    else
        delta_phi_deg = 90;
    end;
    amplitude_cluster = 1;
    phi_deg = input(['Mean Angle of Incidence at Node B (in' ...
		    ' degrees)? ']);
    AS_deg = -1;
    while (AS_deg < 0)
      AS_deg = input('Azimuth Spread (in degrees)? ');
    end;
  else
    delta_phi_deg     = zeros(1,cluster_number);
    amplitude_cluster = zeros(1,cluster_number);
    phi_deg           = zeros(1,cluster_number);
    AS_deg            = zeros(1,cluster_number);
    string = sprintf('\n');
    disp(string);
    for k = 1:1:cluster_number
        if ~(PAS_type == 1)
            delta_phi_deg(k) = -1;
            while (delta_phi_deg(k) < 0)
                delta_phi_deg(k) = input(['Cluster ',num2str(k),'/', ...
                        num2str(cluster_number),' : Half' ...
                        ' Domain Definition (in degrees)? ']);
            end;
        else
            delta_phi_deg(k) = 90;
        end;
        amplitude_cluster(k) = -1;
        while (amplitude_cluster(k) < 0)
            amplitude_cluster(k) = input(['Cluster ',num2str(k),'/', ...
		    num2str(cluster_number),' : Linear' ...
		    ' Power? ']);
        end;
        phi_deg(k) = input(['Cluster ',num2str(k),'/', ...
		    num2str(cluster_number),' : mean Angle' ...
		    ' of Incidence at Node B (in degrees)? ']);
        AS_deg(k) = -1;
        while (AS_deg(k) < 0)
        	AS_deg(k) = input(['Cluster ',num2str(k),'/', ...
		    num2str(cluster_number),' : Azimuth' ...
		    ' Spread (in degrees)? ']);
        end;
        string = sprintf('\n');
        disp(string);
    end;
  end;
end;