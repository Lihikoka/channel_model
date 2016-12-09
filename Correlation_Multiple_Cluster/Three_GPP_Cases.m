function [PDP_linear, RTx, RRx, Fading_Type, Rice_matrix, DirectionOfMovement_rad, AoA_UE_rad, sigma_UE_rad] = ... 
    Three_GPP_Cases(Connection, ID, carrier_frequency_Hz, ... 
                    NumberOfAntennas_NodeB, Spacing_NodeB, ... 
                    NumberOfAntennas_UE, Spacing_UE) 
 
% function [PDP_linear, RTx, RRx, Fading_Type, Rice_matrix, ... 
%           DirectionOfMovement_rad, AoA_UE_rad, sigma_UE_rad] = ... 
%           Three_GPP_Cases(Connection, ID, carrier_frequency_Hz, ... 
%                           NumberOfAntennas_NodeB, Spacing_NodeB, ... 
%                           NumberOfAntennas_UE, Spacing_UE) 
% 
% Provides the PDP, AoA and DoT information of the chosen 3GPP 
% case, plus the spatial correlation properties and the Rice 
% steering matrix. The I/O variables are listed here after. 
% 
% Inputs 
% 
% * Variable Connection, defines whether the simulated channel is 
%   downlink or uplink 
% * Variable ID, defines the 3GPP Case to be simulated 
% * Variable carrier_frequency_Hz 
% * Variable NumberOfAntennas_NodeB, number of antenna elements 
%   of the ULA at Node B 
% * Variable Spacing_NodeB, spacing between the antenna elements 
%   of the ULA at Node B 
% * Variable NumberOfAntennas_UE, number of antenna elements 
%   of the ULA at UE 
% * Variable Spacing_UE, spacing between the antenna elements 
%   of the ULA at UE 
% 
% Outputs 
% 
% * 2-D matrix PDP_linear of size 2 x NumberOfPaths, containing 
%   the linear values of the taps power gains, normalised to 1 
%   and the absolute delays of the taps in seconds 
% * 3-D matrix RTx of size NumberOfPaths x NumberOfTxAntennas 
%   x NumberOfTxAntennas, containing the correlation coefficients 
%   at the ULA transmit antenna with SpacingTx 
% * 3-D matrix RRx of size NumberOfPaths x NumberOfRxAntennas 
%   x NumberOfRxAntennas, containing the correlation coefficients 
%   at the ULA receive antenna with SpacingRx 
% * String Fading_type indicating the nature of the Doppler spectrum 
% * 3-D matrix Rice_amtrix of size NumberOfPaths x NumberOfRxAntennas 
%   x NumberOfTxAntennas whose elements are the product of the phasors 
%   related to the elements of two ULA (Uniform Linear Array), one at 
%   Tx and one at Rx, under incidence angles AoA_Tx_rad and AoA_Rx_rad 
%   respectively. Each of the NumberOfPaths 2-D submatrices is weighted 
%   by the associated K factor of the path 
% * Vector DirectionOfMovement_UE_rad of size 1 x NumberOfPaths 
%   giving the direction of travel of the UE relative to UE antenna 
%   orientation. 
% * Vector AoA_UE_rad of size 1 x NumberOfPaths giving the angle of 
%   arrival of the waves impinging at UE relative to UE antenna 
%   orientation. 
% * Vector sigma_UE_rad of size 1 x NumberOfPaths giving the standard 
%   deviation of the waves impinging at UE relative to UE antenna 
%   orientation. 
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
 
switch ID 
case 1 
    % Rayleigh distributed 
    PDP_dB = [0;  % Average power [dB] 
              0]; % Relative delay (ns) 
           
    % Node B 
    AoA_NodeB_deg = 0; 
    AS_NodeB_deg  = -1; % Irrelevant parameter in Case 1 
    Type_NodeB    = -1; % Irrelevant parameter in Case 1 
    % UE 
    DirectionOfMovement_deg = 0; 
    AoA_UE_deg              = 0; 
    AS_UE_deg               = -1; % Irrelevant parameter in Case 1 
    Type_UE                 = -1; % Irrelevant parameter in Case 1 
    % Fading 
    Fading_Type = 'classic'; 
    % Rice 
    K_factor_dB = -100; 
case 2 
    % ITU Pedestrian A 
    PDP_dB = [0  -9.7  -19.2  -22.8;   % Average power [dB] 
              0 110e-9 190e-9 410e-9]; % Relative delay (ns) 
    % Node B 
    AoA_NodeB_deg = 20.*ones(1,4); 
    AS_NodeB_deg  = 5.*ones(1,4); 
    Type_NodeB    = 3.*ones(1,4); 
    % UE 
    DirectionOfMovement_deg = zeros(1,4); 
    AoA_UE_deg              = 22.5.*ones(1,4); 
    AS_UE_deg               = (180/sqrt(3)).*ones(1,4); 
    Type_UE                 = ones(1,4); 
    % Fading 
    Fading_Type = 'classic'; 
    % Rice 
    K_factor_dB = [3, -100, -100, -100]; 
case 3 
    % ITU Vehicular A 
    PDP_dB = [0  -1     -9     -10     -15     -20;     % Average power [dB] 
              0 310e-9 710e-9 1090e-9 1730e-9 2510e-9]; % Relative delay (ns) 
    % Node B 
    AoA_NodeB_deg = 20.*ones(1,6); 
    AS_NodeB_deg  = 10.*ones(1,6); 
    Type_NodeB    = 3.*ones(1,6); 
    % UE 
    DirectionOfMovement_deg = 22.5.*ones(1,6); 
    AoA_UE_deg              = 67.5.*ones(1,6); 
    AS_UE_deg               = 35.*ones(1,6); 
    Type_UE                 = 3.*ones(1,6); 
    % Fading 
    Fading_Type = 'laplacian'; 
    % Rice 
    K_factor_dB = [-100, -100, -100, -100, -100, -100]; 
case 4 
    % ITU Pedestrian B 
    PDP_dB = [0  -.9   -4.9    -8.0    -7.8   -23.9;    % Average power [dB] 
              0 200e-9 800e-9 1200e-9 2300e-9 3700e-9]; % Relative delay (ns) 
    % Node B 
    AoA_NodeB_deg = [2 -20 10 -8 -33 31]; 
    AS_NodeB_deg  = 15.*ones(1,6); 
    Type_NodeB    = 3.*ones(1,6); 
    % UE 
    DirectionOfMovement_deg = -22.5.*ones(1,6); 
    AoA_UE_deg              = [22.5, -67.5, 22.5, -67.5, 22.5, -67.5]; 
    AS_UE_deg               = 35.*ones(1,6); 
    Type_UE                 = 3.*ones(1,6); 
    % Fading 
    Fading_Type = 'laplacian'; 
    % Rice 
    K_factor_dB = [-100, -100, -100, -100, -100, -100]; 
otherwise 
    disp('Undefined case. Exiting...'); 
    return; 
end; 
% PDP in linear values 
PDP_linear = [10.^(.1.*PDP_dB(1,:)); 
              PDP_dB(2,:)]; 
% Normalisation of PDP 
PDP_linear(1,:) = PDP_linear(1,:)./sum(PDP_linear(1,:)); 
% Direction of movement 
DirectionOfMovement_rad = DirectionOfMovement_deg*pi/180; 
% Impact of connection direction 
if (strcmp(Connection, 'downlink')) 
    % Tx 
    AoA_Tx_deg         = AoA_NodeB_deg; 
    AS_Tx_deg          = AS_NodeB_deg; 
    Type_Tx            = Type_NodeB; 
    NumberOfTxAntennas = NumberOfAntennas_NodeB; 
    SpacingTx          = Spacing_NodeB; 
    % Rx 
    AoA_Rx_deg         = AoA_UE_deg; 
    AoA_UE_rad         = AoA_UE_deg.*pi./180; 
    AS_Rx_deg          = AS_UE_deg; 
    Type_Rx            = Type_UE; 
    NumberOfRxAntennas = NumberOfAntennas_UE; 
    SpacingRx          = Spacing_UE; 
else 
    % Tx 
    AoA_Tx_deg         = AoA_UE_deg; 
    AoA_UE_rad         = AoA_UE_deg.*pi./180; 
    AS_Tx_deg          = AS_UE_deg; 
    Type_Tx            = Type_UE; 
    NumberOfTxAntennas = NumberOfAntennas_UE; 
    SpacingTx          = Spacing_UE; 
    % Rx 
    AoA_Rx_deg         = AoA_NodeB_deg; 
    AS_Rx_deg          = AS_NodeB_deg; 
    Type_Rx            = Type_NodeB; 
    NumberOfRxAntennas = NumberOfAntennas_NodeB; 
    SpacingRx          = Spacing_NodeB; 
end; 
% Derivation of correlation 
if (ID ==1) 
    RTx(1,:,:)   = eye(NumberOfTxAntennas); 
    RRx(1,:,:)   = eye(NumberOfRxAntennas); 
    sigma_UE_rad = -1; % Irrelevant parameter in Case 1 
else 
    for (ii = 1:size(PDP_dB,2)) 
        [temp_RTx, QTx, sTx_deg] = correlation(NumberOfTxAntennas, SpacingTx, ... 
            linspace(0,(NumberOfTxAntennas-1)*SpacingTx, NumberOfTxAntennas), 1, 1,... 
            Type_Tx(ii), AoA_Tx_deg(ii), AS_Tx_deg(ii), 180, 0); 
        [temp_RRx, QRx, sRx_deg] = correlation(NumberOfRxAntennas, SpacingRx, ... 
            linspace(0,(NumberOfRxAntennas-1)*SpacingRx, NumberOfRxAntennas), 1, 1,... 
            Type_Rx(ii), AoA_Rx_deg(ii), AS_Rx_deg(ii), 180, 0); 
        RTx(ii,:,:)     = temp_RTx; 
        RRx(ii,:,:)     = temp_RRx; 
        if strcmp(Connection,'downlink') 
            sigma_UE_rad(1,ii) = sRx_deg * pi/180; 
        else 
            sigma_UE_rad(1,ii) = sTx_deg * pi/180; 
        end; 
    end; 
end; 
% Initialisation of the Rice matrix 
Rice_matrix = init_Rice(K_factor_dB, carrier_frequency_Hz, ... 
    NumberOfTxAntennas, SpacingTx, AoA_Tx_deg*pi/180, ... 
    NumberOfRxAntennas, SpacingRx, AoA_Rx_deg*pi/180);
