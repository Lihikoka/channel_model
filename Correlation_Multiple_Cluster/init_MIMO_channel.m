function [C, R] = init_MIMO_channel(RTx,RRx,Type) 
 
% function C = init_MIMO_channel(RTx,RRx,Type) 
% 
% Performs the Cholesky (Type = 'complex') or the Square-Root Matrix 
% Decomposition (Type = 'real') of the Kronecker product of 
% the (RTx, RRx) pair associated with each path. 
%  
% Inputs 
% 
% * 3-D matrix RTx of size NumberOfPaths x NumberOfTxAntennas 
%   x NumberOfTxAntennas whose elements represent the spatial 
%   correlation properties at Tx 
% * 3-D matrix RRx of size NumberOfPaths x NumberOfRxAntennas 
%   x NumberOfRxAntennas whose elements represent the spatial 
%   correlation properties at Rx 
% * Variable Type, defines the nature of the spatial correlation 
%   properties (complex/real) and the subsequent decomposition 
%   (Cholesky/Square-Root Matrix) 
% 
% Outputs 
%  
% * 2-D block diagonal matrix C of size (NumberOfPaths 
%   * NumberOfTxAntennas * NumberOfTxAntennas) x (NumberOfPaths 
%   * NumberOfRxAntennas * NumberOfRxAntennas), containing 
%   NumberOfPaths spatial correlation shaping matrices of size 
%   (NumberOfTxAntennas * NumberOfRxAntennas) x (NumberOfTxAntennas 
%   * NumberOfRxAntennas) produced by the chosen decomposition 
% * 2-D block diagonal matrix R of size (NumberOfPaths 
%   * NumberOfTxAntennas * NumberOfTxAntennas) x (NumberOfPaths 
%   * NumberOfRxAntennas * NumberOfRxAntennas), containing 
%   NumberOfPaths spatial correlation matrices of size 
%   (NumberOfTxAntennas * NumberOfRxAntennas) x (NumberOfTxAntennas 
%   * NumberOfRxAntennas) 
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
% including, but not limited to, direct, indirect, special, or 
% consequential damages arising out of, resulting from, or any way 
% connected to the use of the item, whether or not based upon 
% warranty, contract, tort, or otherwise; whether or not injury was 
% sustained by persons or property or otherwise; and whether or not 
% loss was sustained from, or arose out of, the results of, the 
% item, or any services that may be provided by CSys. 
% 
% (c) Laurent Schumacher, AAU-TKN/IES/KOM/CPK/CSys - February 2003 
 
C = zeros(size(RTx,2)*size(RRx,2), size(RTx,2)*size(RRx,2)); 
R = zeros(size(RTx,2)*size(RRx,2), size(RTx,2)*size(RRx,2)); 
for ii = 1:size(RTx,1) 
    span = (ii-1)*size(RTx,2)*size(RRx,2)+1:1:ii*size(RTx,2)*size(RRx,2); 
    temp_RTx = reshape(RTx(ii,:,:),size(RTx,2),size(RTx,3)); 
    temp_RRx = reshape(RRx(ii,:,:),size(RRx,2),size(RRx,3)); 
    if (Type == 'complex') 
        temp_C = chol(kron(temp_RTx, temp_RRx))'; 
    else 
        temp_C = sqrtm(sqrt(kron(temp_RTx, temp_RRx))); 
        if (max(max(abs(imag(temp_C)))) > 0) 
            disp('Correlation matrix non-positive definite. Exiting...'); 
            break; 
        end; 
    end; 
    R(span, span) = kron(temp_RTx, temp_RRx); 
    C(span, span) = temp_C; 
end;