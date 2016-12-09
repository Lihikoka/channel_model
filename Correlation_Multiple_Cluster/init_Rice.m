function y = Rice_matrix(K_factor_dB, carrier_frequency_Hz, ... 
    NumberOfTxAntennas, spacing_Tx, AoD_Tx_rad, ... 
    NumberOfRxAntennas, spacing_Rx, AoA_Rx_rad) 
 
% function y = Rice_matrix(K_factor_dB, carrier_frequency_Hz, ... 
%                          NumberOfTxAntennas, spacing_Tx, AoD_Tx_rad, ... 
%                          NumberOfRxAntennas, spacing_Rx, AoA_Rx_rad) 
% 
% Computes the Rice steering matrix as described in Appendix B of 
% 3GPP document R1-01-1179. 
% 
% Inputs 
% 
% * Variable K_factor_dB, Rice factor in dB 
% * Variable carrier_frequency_Hz, carrier frequency in Hz 
% * Variable NumberOfTxAntennas, number of antenna elements of 
%   the Uniform Linear Array (ULA) at Tx 
% * Variable spacing_Tx, spacing of the antenna elements of 
%   the ULA at Tx 
% * Variable AoD_Tx_rad, angle of departure of the waves at Tx 
% * Variable NumberOfRxAntennas, number of antenna elements of 
%   at Rx 
% * Variable spacing_Rx, sapcing of the antenna elements of 
%   the ULA at Rx 
% * Variable AoA_Rx_rad, angle of arrival of the waves at Rx 
% 
% Output 
% 
% A 3-D matrix of size NumberOfPaths x NumberOfRxAntennas 
% x NumberOfTxAntennas whose elements are the product of the phasors 
% related to the elements of two ULAs, one at Tx and one at Rx, under 
% incidence angles AoD_Tx_rad and AoA_Rx_rad respectively. Each of 
% the NumberOfPaths 2-D submatrices is weighted by the associated 
% K factor of the path. 
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
% (c) Laurent Schumacher, AAU-TKN/IES/KOM/CPK/CSys - February 2002 
 
wavelength_m = 3e8/carrier_frequency_Hz; 
y            = zeros(size(K_factor_dB,2), NumberOfRxAntennas, NumberOfTxAntennas); 
 
for ii=1:size(K_factor_dB,2) 
    K            = 10^(.1*K_factor_dB(ii)); 
    step_Tx      = exp(j*2*pi*spacing_Tx*sin(AoD_Tx_rad(ii))/wavelength_m); 
    vector_Tx    = step_Tx.^(0:NumberOfTxAntennas-1); 
    step_Rx      = exp(j*2*pi*spacing_Rx*sin(AoA_Rx_rad(ii))/wavelength_m); 
    vector_Rx    = step_Rx.^((0:NumberOfRxAntennas-1).'); 
    y(ii,:,:)    = sqrt(K).*(vector_Rx * vector_Tx); 
end;