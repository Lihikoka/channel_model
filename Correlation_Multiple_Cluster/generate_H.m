function H = generate_H(nTx, nRx, R, lH, type, PDP_linear) 
% H = generate_H(nTx, nRX, R, lH, type) 
% 
% Generates a H matrix whose correlation properties fit the values 
% in given in correlation matrix R. These values can be (complex) 
% field correlation coefficients or (real positive) power 
% correlation coefficients 
 
if (type == 0) % complex correlation 
    C = sqrtm(R);%chol(R); 
else % power correlation coefficient 
    C  = sqrtm(R.^.5); 
end; 
 
H  = zeros(nRx, nTx, lH); 
for k = 1:lH 
  TildeA   = (C*(randn(nTx*nRx,1)+i*randn(nRx*nTx,1)))./sqrt(2); 
  H(:,:,k) = reshape(TildeA, nRx, nTx); 
end;