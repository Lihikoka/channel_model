function clusterPowerRatios = getClusterPowers(tau_n,Gamma,sigmaCluster)
% Generate time cluster powers (relative to 1 mW).
%
% Inputs:
%   - tau_n: array of time cluster excess delays, in ns
%   - Gamma: time cluster decay constant, in ns
%   - sigmaCluster: per-cluster shadowing, in dB
% Output:
%   - clusterPowers: array of time cluster powers, relative to 1 mW
%
% Copyright © 2016 NYU

%%% number of clusters
numberOfTimeClusters = size(tau_n,2);

%%% generate per-cluster shadowing
Z = sigmaCluster*randn([1 numberOfTimeClusters]); 

%%% cluster ratios
clusterPowerRatios_temp = exp(-tau_n/Gamma).*10.^(Z/10);

%%% normalize cluster ratios such that their sum equals 1
clusterPowerRatios = clusterPowerRatios_temp/sum(clusterPowerRatios_temp);

end