function [Q, sigma_deg] = normalisation_laplacian(number_clusters, ...
						  power_lin, AS_deg, delta_phi_deg)
% [Q, sigma_deg] = normalisation_laplacian(number_clusters,
%                  power_lin, AS_deg, delta_phi_deg)
%
% Computes the variances sigma_k and the power normalising
% coefficients Q_k such that the Power Azimuth Spectrum (PAS) can
% be regarded as a probability distribution function (pdf), that is
% to say int_{-\pi}^{pi} PAS(\phi) d\phi = 1.
%
% Revision history
%
% July 2003  - Bug fix: half-domain definition delta_phi_deg held as
%              a vector instead of as a scalar
% April 2003 - Bug fix: half-domain definition delta_phi_deg
%              passed as parameter. Used to be hard-coded. Thanks to
%              Cesar Caballero, University of Zaragoza (Spain) for
%              reporting the bug.
% July 2001  - Creation
%
%
% STANDARD DISCLAIMER
%
% The Computer Science Institute of the University of Namur (hereafter
% "FUNDP-INFO") is furnishing this item "as is". FUNDP-INFO does not
% provide any warranty of the item whatsoever, whether express,
% implied, or statutory, including, but not limited to, any warranty
% of merchantability or fitness for a particular purpose or any
% warranty that the contents of the item will be error-free.
%
% In no respect shall FUNDP-INFO incur any liability for any damages,
% including, but limited to, direct, indirect, special, or
% consequential damages arising out of, resulting from, or any way
% connected to the use of the item, whether or not based upon
% warranty, contract, tort, or otherwise; whether or not injury was
% sustained by persons or property or otherwise; and whether or not
% loss was sustained from, or arose out of, the results of, the
% item, or any services that may be provided by FUNDP-INFO.
%
% (c) Laurent Schumacher, FUNDP-INFO - July 2003

%
% Computation of sigma
%
for k=1:number_clusters
    string = ['exist(''AS2sigma_laplacian_', num2str(delta_phi_deg(k)),'.mat'')'];
    if eval(string) == 0
        disp(['!!! Look-up table AS -> sigma unavailable for ', ...
                num2str(delta_phi_deg(k)),'¢X half-domain definition']);
        disp('!!! Please generate with AS_vs_sigma_laplacian.m');
        break;
    end;
    string = ['load AS2sigma_laplacian_', num2str(delta_phi_deg(k)),' tmp;'];
    eval(string);
    if ((AS_deg(k) < min(tmp(:,1))) | (AS_deg(k) >= max(tmp(:,1))))
        disp('!!! Look-up table AS -> sigma incomplete.');
        disp('!!! Please regenerate with AS_vs_sigma_laplacian.m');
        break;
    end;
    pos          = find(tmp(:,1)>=AS_deg(k));
    sigma_deg(k) = ((tmp(pos(1),1)-AS_deg(k))*tmp(pos(1)-1,2) + ...
        (AS_deg(k) - tmp(pos(1)-1,1))*tmp(pos(1),2))/ ...
        (tmp(pos(1),1)-tmp(pos(1)-1,1));
end;
% Conversion degree -> rad
AS_rad    = AS_deg.*(pi/180);
sigma_rad = sigma_deg.*(pi/180);
% Domain definition
delta_phi_rad = delta_phi_deg.*(pi/180);
%
% Computation of Q
%
if (number_clusters == 1)
  Q = 1/(1-exp(((-1)*sqrt(2)*delta_phi_rad)/sigma_rad));
else
  A = zeros(number_clusters);
  A(1:number_clusters-1,1) = 1/(sigma_rad(1)*power_lin(1));
  for k=2:number_clusters
    A(k-1,k) = (-1)/(sigma_rad(k)*power_lin(k));
  end;
  A(number_clusters,:) = ones(1,number_clusters)-exp(((-1).* ...
						 sqrt(2).* ...
						 delta_phi_rad)./ ...
						 sigma_rad);
  b = zeros(number_clusters,1);
  b(number_clusters,1) = 1;
  Q = (inv(A)*b).';
end;
%
% Validation
%
if ~(sum((ones(1,number_clusters)-exp(((-1).*sqrt(2).* ...
				       delta_phi_rad)./sigma_rad))*Q.')-1 ...
     < 1e-15)
  disp('Normalisation of laplacian distribution failed!');
end;