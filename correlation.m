function [R, Q, sigma_deg] = correlation(M, spacing, d_norm, ... 
					 cluster_number, amplitude_cluster, ... 
                     PAS_type, phi_deg, AS_deg, ... 
                     delta_phi_deg, type) 
% [R, Q, sigma_deg] = correlation(M, spacing, d_norm, 
%                     cluster_number, amplitude_cluster, PAS_type, 
%                     phi_deg, AS_deg, delta_phi_deg, type) 
% 
% Derives the correlation matrix R from the description of the 
% environment, namely the antenna element spacing and the way the 
% waves impinge (number of clusters, PAS type, AS and mean angle of 
% incidence). The generated correlation coefficients are either 
% (complex) field correlation coefficients (type = 0) or (real 
% positive) power correlation coefficients (type = 1). 
% 
% Inputs 
% 
% * Variable M, number of antenna elements of the ULA 
% * Variable spacing, spacing of the antenna elements of the ULA 
% * Vector d_norm, containing the relative spacings of the M 
%   elements of the ULA with respect to the first one 
% * Variable cluster_number, number of impinging clusters 
% * Vector amplitude_cluster, containing the amplitude of 
%   the cluster_number impinging clusters 
% * Variable PAS_type, defines the nature of the PAS, isotropic 
%   uniform (1), isotropic Gaussian (2), isotropic Laplacian (3) 
%   or directive Laplacian (6) according to 3GPP-3GPP2 SCM AHG 
%   radiation pattern 
% * Vector phi_deg, containing the AoAs of the number_cluster 
%   clusters 
% * Vector AS_deg, containing the ASs of the number_cluster 
%   clusters 
% * Vector delta_phi_deg, containing the constraints limits 
%   of the truncated PAS, if applicable (Gaussian and Laplacian 
%   case) 
% * Variable type, defines whether correlation properties 
%   should be computed in field (0) or in power (1) 
% 
% Outputs 
% 
% * 2-D Hermitian matrix R of size M x M, whose elements are 
%   the correlation coefficients of the corresponding elements 
%   of the ULA impinged by the described PAS 
% * Vector Q of cluster_number elements, giving the normalisation 
%   coefficients to be applied to the clusters in order to have 
%   the PAS fulfilling the definition of a pdf 
% * Vector sigma_deg of cluster_number of elements, giving 
%   the standard deviations of the clusters, derived from 
%   their description. There is not necessarily equality between 
%   the given AS and the standard deviation of the modelling PAS. 
% 
% Revision history 
% 
% April 2003    - Addition of case 6 (directive Laplacian) 
% February 2002 - Creation 
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
% (c) Laurent Schumacher, AAU-TKN/IES/KOM/CPK/CSys - February 2002 
 
% Normalisation 
switch(PAS_type) 
    case 1 % Uniform distribution 
        Q = normalisation_uniform(cluster_number, ... 
					 amplitude_cluster, ... 
					 AS_deg); 
        sigma_deg = -1; % meaningless in uniform PAS 
    case 2 % Gaussian distribution 
        [Q, sigma_deg] = ... 
	    normalisation_gaussian(cluster_number, ... 
				   amplitude_cluster, ... 
				   AS_deg, ... 
                   delta_phi_deg); 
    case {3, 6} % Laplacian distribution 
        [Q, sigma_deg] = ... 
	    normalisation_laplacian(cluster_number, ... 
				    amplitude_cluster, ... 
				    AS_deg, ... 
                    delta_phi_deg); 
    otherwise 
	disp('PAS type non supported. Exiting...'); 
    
end; 
 
if (M>1) 
    % Complex correlation computation 
    switch(PAS_type) 
    case 1 % Isotropic Uniform distribution 
        Rxx = besselj(0,2.*pi.*d_norm); 
        Rxy = zeros(size(Rxx)); 
        for k=1:cluster_number 
            Rxx = Rxx + Q(k).*Rxx_uniform(d_norm, ... 
						  phi_deg(k), ... 
						  AS_deg(k)); 
            Rxy = Rxy + Q(k).*Rxy_uniform(d_norm, ... 
						  phi_deg(k), ... 
						  AS_deg(k)); 
        end; 
    case 2 % Isotropic Gaussian distribution 
        Rxx = besselj(0,2.*pi.*d_norm); 
        Rxy = zeros(size(Rxx)); 
        for k=1:cluster_number 
            Rxx = Rxx + Q(k).*Rxx_gaussian(d_norm, ... 
						  phi_deg(k), ... 
						  sigma_deg(k), ... 
                          delta_phi_deg(k)); 
            Rxy = Rxy + Q(k).*Rxy_gaussian(d_norm, ... 
						  phi_deg(k), ... 
						  sigma_deg(k), ... 
                          delta_phi_deg(k)); 
        end; 
    case 3 % Isotropic Laplacian distribution 
        Rxx = besselj(0,2.*pi.*d_norm); 
        Rxy = zeros(size(Rxx)); 
        for k=1:cluster_number 
            Rxx = Rxx + Q(k).*Rxx_laplacian(d_norm, ... 
						  phi_deg(k), ... 
						  sigma_deg(k), ... 
                          delta_phi_deg(k)); 
            Rxy = Rxy + Q(k).*Rxy_laplacian(d_norm, ... 
						  phi_deg(k), ... 
						  sigma_deg(k), ... 
                          delta_phi_deg(k)); 
        end; 
    case 6 % Directive Laplacian distribution 
        % Hard-coded radiation pattern characteristics 
        tmp = rho_lpl_SCM_AHG_approx(d_norm, ... 
						  phi_deg, ... 
						  Q, ... 
                          sigma_deg, ... 
                          delta_phi_deg, ... 
                          12, 20, 70); 
        Rxx = real(tmp); 
        Rxy = imag(tmp); 
    end; 
    % Correlation coefficient computation 
    if (type == 0) 
        tmp = Rxx + i.* Rxy; 
    else 
        tmp = Rxx.^2 + Rxy.^2; 
    end; 
    R = toeplitz(tmp); 
else 
    R = 1; 
end;
