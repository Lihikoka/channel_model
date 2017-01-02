function [Pr_dBm, PL] = getRXPower(f_str,n,SF,TXPower,dist,d0,dynamicRange)
% Generate the omnidirectional received power and path loss
%
% Inputs:
%   - f_str: a string specifying the frequency of interest
%   - n: the frequency-dependent path loss exponent
%   - SF: the shadow factor, in dB
%   - TXPower: the transmit power, typically set to 0 dBm
%   - dist: the T-R separation, in meters
%   - d0: the free space reference distance, typically 1 meter
%   - dynamicRange: the maximum allowable path loss, typically 180 dB in
%   the NYU measurements
% Output:
%   - Pr_dBm: the omnidirectional received power, in dBm
%   - PL: the omnidirectional path loss, in dB
%
% Copyright © 2016 NYU

switch f_str
    case '28_GHz'
        f_ = 28e9;
    case '73_GHz'
        f_ = 73e9;
    case '2873_GHz'
        f_ = 28e9; % For the joint scenario, use the 28 GHz frequency and path loss parameters
    otherwise
end

% constants
c = 3e8; %% speed of light (m/s)
lambda = c/f_; %% wavelength (m)

% free space path loss at d0 (dB)
PLref = 20*log10(4*pi*d0/lambda);

% absolute path loss at distance dist
PL = min(PLref + n*10*log10(dist/d0)+SF*randn,dynamicRange);

% total received power (dBm) at distance dist 
Pr_dBm = TXPower - PL;

end