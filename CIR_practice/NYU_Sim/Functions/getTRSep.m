function dist = getTRSep(dmin,dmax)
% Generate the T-R separation distance, in meters
%
% Inputs:
%   - dmin: the minimum distance, in meters
%   ns
%   - dmax: the maximum distance, in meters
%   - gamma: subpath decay constant, in ns
% Output:
%   - dist: the T-R separation distance, in meters
%
% Copyright © 2016 NYU

% Generates a uniform random variable distributed between dmin and dmax
dist = dmin + (dmax-dmin).*rand;

end