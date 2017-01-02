% NYU-Simulator Version 1.1 (NYU-Sim), developed by:
%
% Mathew K. Samimi, Shu Sun - NYU WIRELESS, 16 February 2016
%
% This script runs the NYU-Simulator
%
% Variable inputs to this script:
%   N (number of CIRs to genearate)
%   See loadUserInputs.m in the ReadFiles folder
%
% Outputs of this script:
% ['Results_',frequency,'_',environment,'_',scenario,'_',num2str(N),'_CIRs.mat']
% ['SecondOrderStats_',frequency,'_',environment,'_',scenario,'_',num2str(N),'_CIRs.txt'];
%
% Copyright ?2016 NYU
%
% See REFERENCES below for more information:
%
% [1] M. Samimi et al., "28 GHz Angle of Arrival and Angle of Departure 
% Analysis for Outdoor Cellular Communications Using Steerable Beam 
% Antennas in New York City," in 2013 IEEE  Vehicular Technology 
% Conference (VTC Spring), pp.1-6, 2-5 June 2013.
% URL: http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=6691812&isnumber=6691801
% 
% [2] M. K. Samimi, T. S. Rappaport, "Ultra-wideband statistical channel
% model for non line of sight millimeter-wave urban channels," in 2014
% IEEE Global Communications Conference (GLOBECOM), pp. 3483-3489, 8-12
% Dec. 2014.
% URL: http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=7037347&isnumber=7036769
% 
% [3] M. K. Samimi, T. S. Rappaport, "3-D statistical channel model for
% millimeter-wave outdoor mobile broadband communications," 2015 IEEE 
% International Conference on in Communications (ICC), pp.2430-2436, 
% 8-12 June 2015.
% URL: http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=7248689&isnumber=7248285
% 
% [4] M. K. Samimi, T. S. Rappaport, “Statistical Channel Model with 
% Multi-Frequency and Arbitrary Antenna Beamwidth for Millimeter-Wave
% Outdoor Communications,?in 2015 IEEE Global Communications Conference,
% Exhibition & Industry Forum (GLOBECOM) Workshop, Dec. 6-10, 2015.
% URL: http://arxiv.org/abs/1510.03081
% 
% [5] M. K. Samimi, T. S. Rappaport, ?8 GHz Millimeter-Wave Ultrawideband 
% Small-Scale Fading Models in Wireless Channels,?in 2016 IEEE 
% Vehicular Technology Conference (VTC2016-Spring), 15-18 May, 2016.
% URL: http://arxiv.org/abs/1511.06938
% 
% [6] M. K. Samimi, S. Sun, and T. S. Rappaport, “MIMO Channel Modeling 
% and Capacity Analysis for 5G Millimeter-Wave Wireless Systems,?in 
% the 10th European Conference on Antennas and Propagation (EuCAP?016), 
% April 2016.
% URL: http://arxiv.org/abs/1511.06940
% 
% [7] M. K. Samimi, T. S. Rappaport, “Local Multipath Model Parameters for 
% Generating 5G Millimeter-Wave 3GPP-like Channel Impulse Response,?
% in the 10th European Conference on Antennas and Propagation (EuCAP?016),
% April 2016.
% URL: http://arxiv.org/abs/1511.06941
 
clear 
clc
close all

% master folder, and add path to all underlying subfolders
setMasterFolder

% Folder where the output is saved;
resultsFolder = [masterFolder,'\Results'];

% Load user inputs
loadUserInputs

%%%%%%%%%%%%%%%%%
% NYU-Simulator %
%%%%%%%%%%%%%%%%%

% Generate N Channel Impulse Reponses (CIRs): one for each of the N users
% Output of generateCIR: CIR_Struct
generateCIR

% Output the generated CIRs
CIR_Struct;

% Generate the local area CIRs based on the previously generated single
% CIRs
CIR_Struct = getLocalCIR(CIR_Struct,TxArrayType,RxArrayType,Nt,Nr,Wt,Wr,dTxAnt,dRxAnt);

% Compute RMS DS and RMS Angular Spreads, and output to .txt files
cd(resultsFolder)
extractStatistics(CIR_Struct)

% Save output for later processing
save(['Results_',frequency,'_',scenario,'_',num2str(N),'_CIRs'],'CIR_Struct')

%===
% Plotting
h1_delay = CIR_Struct.CIR_1.pathDelays;
h1_power = CIR_Struct.CIR_1.pathPowers;

h2_delay = CIR_Struct.CIR_2.pathDelays;
h2_power = CIR_Struct.CIR_2.pathPowers;

subplot(1,2,1);
stem(h1_delay, h1_power, 'r');
xlabel('Path delay (ns)');
ylabel('Path power (mW)');

subplot(1,2,2);
stem(h2_delay, h2_power, 'b');
xlabel('Path delay (ns)');
ylabel('Path power (mW)');















