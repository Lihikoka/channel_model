clear all; 
close all; 
 
% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
% The ones need to be changed in TLT-6206 project work % 
% Antenna configuration at Node B (for the capacity analysis to work, keep nTx and nRx equal) 
NumberOfAntennas_NodeB = 2; 
Spacing_NodeB          = 0.5; % antenna spacing in wavelengths 
% Antenna configuration at UE (for the capacity analysis to work, keep nTx and nRx equal) 
NumberOfAntennas_UE = 2; 
Spacing_UE          = 0.5; % antenna spacing in wavelengths 
% 3GPP case to be simulated 
Three_GPP_Case = 2; % Channel model, see Three_GPP_Cases.m 
% Number of runs of the main loop (256, 512, 1024, ...) increases statistical
% reliability but also the simulation time. Use at least 256. 
NumberOfIterations = 256 ; 
 
% The parameters from here you don't have to touch % 
% Direction of connection 
Connection = 'downlink'; % let's stay in DL Direction 
% Chip rate 
ChipRate_Hz = 3.84e6; 
% Carrier frequency 
CarrierFrequency_Hz = 2.15e9; 
% Speed 
Speed_kmh = 120; 
CorrelationCoefficientType = 'complex'; % Complex field or Real power correlation coefficients 
% NumberOfIterations of the vector of fading coefficients 
FadingNumberOfIterations = 1024; 
% Oversampling factor of the Doppler spectrum 
FadingOversamplingFactor = 4; 
% Number of chips dealt with during one iteration 
NumberOfChipsPerIteration = 512;  
% Chip oversampling factor 
ChipOversamplingFactor = 1; 
 
 
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
% consequential damages arising out of, resulting from, or any wayconnected to the use of the item, whether or not based upon 
% warranty, contract, tort, or otherwise; whether or not injury was 
% sustained by persons or property or otherwise; and whether or not 
% loss was sustained from, or arose out of, the results of, the 
% item, or any services that may be provided by CSys. 
% 
% (c) Laurent Schumacher, AAU-TKN/IES/KOM/CPK/CSys - February 2002 
 
% Initialisation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
% Parameters' conversion 
 
Speed_ms          = Speed_kmh*3.6; 
Wavelength_m      = 3e8/CarrierFrequency_Hz; 
Max_Doppler_shift = Speed_ms/Wavelength_m; 
 
% Parameters assessment 
 
FadingSamplingTime = Wavelength_m / (2 * FadingOversamplingFactor * Speed_ms); 
FadingLength       = FadingNumberOfIterations * FadingSamplingTime * Speed_ms; 
if (FadingLength < (40*Wavelength_m)) 
    disp('Warning! Length of fading vector smaller than 40 wavelengths. Risk of correlation'); 
end; 
 
% Derivation of the PDP and the correlation matrices 
[PDP_linear, RTx, RRx, FadingType, Rice_matrix, DirectionOfMovement_UE_rad, AoA_UE_rad, sigma_UE_rad] = ... 
    Three_GPP_Cases(Connection, Three_GPP_Case, CarrierFrequency_Hz,... 
                    NumberOfAntennas_NodeB, Spacing_NodeB, ... 
                    NumberOfAntennas_UE, Spacing_UE); 
% Initialisation of size variables 
NumberOfTxAntennas     = size(RTx, 2); 
NumberOfRxAntennas     = size(RRx, 2); 
NumberOfPaths          = size(PDP_linear,2); 
% Computation of the matrix of fading coefficients 
FadingMatrix = init_fading(FadingNumberOfIterations, FadingOversamplingFactor, ... 
                           NumberOfTxAntennas, NumberOfRxAntennas, ... 
			               NumberOfPaths, FadingType, sigma_UE_rad, ... 
                           DirectionOfMovement_UE_rad-AoA_UE_rad); 
% Computation of the correlation matrices, one for each tap 
[C, R] = init_MIMO_channel(RTx, RRx, CorrelationCoefficientType); 
% Spatial correlation of the fading matrix 
FadingMatrix = C * FadingMatrix; 
% Initalisation of normalisation diagonal matrix 
pdp_coef = []; 
for ii = 1:size(PDP_linear,2) 
    pdp_coef = [pdp_coef, sqrt(PDP_linear(1,ii)).*ones(1,size(FadingMatrix,1)/size(PDP_linear,2))]; 
end; 
% Normalisation of the correlated fading processes 
FadingMatrix = diag(pdp_coef)*FadingMatrix; 
% Addition of the Rice component on the first path (3GPP case 2) 
if (Three_GPP_Case == 2) 
    time_vector = exp((j*pi*cos(DirectionOfMovement_UE_rad(1)-AoA_UE_rad(1))/FadingOversamplingFactor) ... 
        .*(0:FadingNumberOfIterations-1)); 
    Rice_vector  = reshape(Rice_matrix(1,:,:), NumberOfTxAntennas*NumberOfRxAntennas, 1); 
    FadingMatrix = FadingMatrix + ... 
        [Rice_vector; zeros((NumberOfPaths-1)*NumberOfTxAntennas*NumberOfRxAntennas,1)]*time_vector; 
    % Sum over all paths equal to 1 
    FadingMatrix = FadingMatrix./sqrt((sum(sum(abs(FadingMatrix).^2))... 
        /NumberOfTxAntennas/NumberOfRxAntennas/size(FadingMatrix,2))); 
    R(1:NumberOfTxAntennas*NumberOfRxAntennas,1:NumberOfTxAntennas*NumberOfRxAntennas) = ... 
        ((Rice_vector*Rice_vector') + R(1:NumberOfTxAntennas*NumberOfRxAntennas,1:NumberOfTxAntennas*NumberOfRxAntennas))... 
        ./(10^.3+1); 
end; 
% Set-up of the iteration process 
NumberOfTransmittedSamples = NumberOfChipsPerIteration; 
NumberOfReceivedSamples    = NumberOfTransmittedSamples + ... 
                             floor(PDP_linear(2,NumberOfPaths)*ChipRate_Hz*ChipOversamplingFactor) ... 
                             +1; 
                          
% Main Loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
RunningTimeInstant = 0; 
% Random choice of the sample to be tracked for the valdiation process 
TrackedSample = 1 + floor((NumberOfChipsPerIteration-1)*rand(1)); 
for ii = 1:NumberOfIterations; 
    time = clock; 
    disp([num2str(time(4)),':',num2str(time(5)),'.',num2str(floor(time(6))),... 
            ' - Iteration ',num2str(ii),'/',num2str(NumberOfIterations)]); 
    [Channel, CorrelatedFading] = MIMO_channel(FadingMatrix, FadingOversamplingFactor, ... 
        Speed_ms, CarrierFrequency_Hz, RunningTimeInstant, ChipOversamplingFactor, ... 
        ChipRate_Hz, NumberOfChipsPerIteration, PDP_linear, NumberOfTxAntennas, ... 
        NumberOfTransmittedSamples, NumberOfRxAntennas, NumberOfReceivedSamples); 
    H(:,ii)            = CorrelatedFading(:,TrackedSample); 
    RunningTimeInstant = RunningTimeInstant + (NumberOfChipsPerIteration/ChipRate_Hz); 
end; 
temp_H = H; 
clear H; 
for ii = 1:NumberOfIterations 
    for jj = 1:NumberOfPaths 
        for kk = 1:NumberOfRxAntennas 
            for ll = 1:NumberOfTxAntennas 
                 H(ll,kk,jj,ii) = temp_H(((jj-1)*NumberOfTxAntennas+(ll-1))*NumberOfRxAntennas+kk, ii); 
            end; 
        end; 
    end; 
end; 
save H; 
 
% Validation of MIMO channel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
plot_MIMO;
