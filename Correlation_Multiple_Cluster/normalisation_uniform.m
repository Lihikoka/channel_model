function Q = normalisation_uniform(number_clusters, power_lin, AS_deg)
% Q = normalisation_uniform(number_clusters, power_lin, AS_deg)
%
% Computes the power normalising coefficients Q_k such that the
% Power Azimuth Spectrum (PAS) can be regarded as a probability
% distribution function (pdf), that is to say int_{-\pi}^{pi}
% PAS(\phi) d\phi = 1.
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

%
% Computation
%
delta_phi_rad = (AS_deg.*pi./180).*sqrt(3);
if (number_clusters == 1)
  Q = 1/(2*delta_phi_rad);
else
  A = zeros(number_clusters);
  A(1:number_clusters-1,1) = 1/power_lin(1);
  for k=2:number_clusters
    A(k-1,k) = (-1)/power_lin(k);
  end;
  A(number_clusters,:) = delta_phi_rad;
  b = zeros(number_clusters,1);
  b(number_clusters,1) = .5;
  Q = (inv(A)*b).';
end;
%
% Validation
%
if ((sum(delta_phi_rad*Q.') - .5) > 1e-9)
  disp('Normalisation of uniform distribution failed!');
end;
