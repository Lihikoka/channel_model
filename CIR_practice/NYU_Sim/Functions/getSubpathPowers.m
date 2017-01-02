function subpathPowers = ...
    getSubpathPowers(rho_mn,clusterPowers,gamma,sigmaSubpath,TXPower,dynamicRange)
% Generate the subpath powers
%
% Inputs:
%   - rho_mn: a structure containing the intra-cluster subpath delays, in
%   ns
%   - clusterPowers: an array containing the time cluster powers, relative
%   to 1 mW
%   - gamma: subpath decay constant, in ns
%   - sigmaSubpath: per-subpath shadowing, in dB
%   - dynamicRange: maximum allowable path loss for each path, in dB
% Output:
%   - subpathPowers: a structure containing the subpath powers, relative to
%   1 mW
%
% Copyright © 2016 NYU

%%% number of clusters
numberOfClusters = size(clusterPowers,2);

%%% initialize the structure that will contain component powers
subpathPowers = struct;   

for clusterIndex = 1:numberOfClusters
    
    %%% current intra-cluster delays
    rho = rho_mn.(['c',num2str(clusterIndex)]);
    
    %%% number of components in current cluster
    numberOfComponents = numel(rho);
    
    %%% per sub path shadowing
    U = sigmaSubpath*randn([1 numberOfComponents]);        
        
    %%% generate sub path ratios
    subPathRatios_temp = exp(-rho/gamma).*10.^(U/10);
    
    %%% cluster power
    clusterPower = clusterPowers(clusterIndex);

    %%% normalize subpath power ratios such that their sum equals 1
    subPathRatios = subPathRatios_temp/sum(subPathRatios_temp)*clusterPower;    
        
    %%% ensure the subpath powers are at least equal to the dynamicRange
    powerTemp = max(subPathRatios,10^((TXPower-dynamicRange)/10));
    
    %%% store sub path powers    
    subpathPowers.(['c',num2str(clusterIndex)]) = powerTemp;
  
    
end



end
