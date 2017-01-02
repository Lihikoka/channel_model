% Copyright ?2016 NYU

% Frequency of interest: 
%   - '2873_GHz' (Joint 28 - 73 GHz statistics)
%   - '28_GHz'
%   - '73_GHz'
frequency = '2873_GHz';

% Scenario of interest: 'Outdoor'
scenario = 'Outdoor';

% RX access type of interest: 'Mobile'
accessType = 'Mobile';

% Number of desired independent users in the current simulation.
N = 20;

% Percentage of LOS users (a number between 0 and 1)
% PercentLOS only works for the Joint scenario, not the individual
% frequency scenarios.
if strcmp(frequency,'2873_GHz') == false
    PercentLOS = 0;
else
    PercentLOS = 1;
end


% Transmit (TX) power in dBm
TXPower = 0;

% Maximum possible path loss (dB)
dynamicRange = 180;

% Antenna gain specifications - antenna half-power beamwidths (degrees)
% Azimuth (\theta) and Elevation (\phi) half-power beamwidth at the 
% transmitter side
theta_3dB_TX = 10; phi_3dB_TX = 10;

% Azimuth (\theta) and Elevation (\phi) half-power beamwidth at the 
% receiver side
theta_3dB_RX = 10; phi_3dB_RX = 10;

% Input parameters for the function getLocalCIR
TxArrayType = 'ULA'; RxArrayType = 'ULA'; 
Nt = 81; Nr = 1; 
dTxAnt = 1/2; dRxAnt = 1/2; 
Wt = 1; Wr = 2;

% Velocity magnitude (m/s) and directions of travel (deg) for each of the N
% users
v = 2*rand(N,1);
theta_v = round(360*rand(N,1));

% Durations of travel (s) for each of the N users. Delta_t = 0 s
% corresponds to one static omnidirectional TX-RX link. A non-zero Delta_t 
% corresponds to a moving user.
Delta_t = randi(10,[N 1]);

% Spatial sampling
dr = .5;


