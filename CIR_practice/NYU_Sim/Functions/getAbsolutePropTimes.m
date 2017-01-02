function t_mn = getAbsolutePropTimes(dist,tau_n,rho_mn)
% Generates absolute time of multipath components.
%
% Inputs:
%   - dist: T-R separation distance in meters
%   - tau_n: array of time cluster delays in ns
%   - rho_mn: structure of intra-cluster subpath delays in ns
% Ouput:
%   - t_mn: structre of intra-cluster subpath absolute times of arrival in
%   ns
%
% Copyright © 2016 NYU

    %%% initialize structure that will contain propagatin times
    t_mn = struct;

    %%% number of clusters
    numberOfClusters = size(tau_n,2);

    %%% speed of light (m/s)
    c = 3e8;
    
    %%% absolute propagation time of first arrival in ns
    t0 = dist/c*1e9;

    for clusterIndex = 1:numberOfClusters

        %%% cluster excess delay       
        tau = tau_n(clusterIndex);

        %%% intra cluster excess delays
        rho = rho_mn.(['c',num2str(clusterIndex)]);

        %%% recover absolute propagation times of arrival
        t_mn.(['c',num2str(clusterIndex)]) = t0+tau+rho;

    end



end