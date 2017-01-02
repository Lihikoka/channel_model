function powerSpectrum = getPowerSpectrum(numberOfClusterSubPaths,t_mn,subPathPowers,phases_mn,...
    subpath_AODs,subpath_AOAs,Pr_dBm)
% Generate the multipath parameters 
%
% Inputs:
%   - numberOfClusterSubPaths: an array containing the number of subpaths
%   in each time cluster
%   - t_mn: a structure containing the time of arrivals of each subpath, in
%   ns
%   - subPathPowers: a structure containing the power levels of the
%   multipaths, relative to 1 mW
%   - phases_mn: a structure containing the phases of the multipaths, in
%   radians
%   - subpath_AODs: a structure containing the (AODs,ZODs) of the
%   multipaths, in degrees
%   - subpath_AOAs: a structure containing the (AOAs,ZOAs) of the
%   multipaths, in degrees 
%   - Pr_dBm: total omnidirectional received power, in dBm
% Output:
%   - powerSpectrum: an array containing the multipath parameters:
%       Column 1: path delays
%       Column 2: path powers, relative to 1 mW
%       Column 3: path phases, in radians
%       Column 4: azimuth angles of departure, in degrees
%       Column 5: zenith angles of departure, in degrees
%       Column 6: azimuth angles of arrival, in degrees
%       Column 6: zenith angles of departure, in degreeshe TX
%
% Copyright © 2016 NYU

%%% number of time clusters
numberOfClusters = numel(numberOfClusterSubPaths);

%%% initialize power spectrum
powerSpectrum = zeros(sum(numberOfClusterSubPaths),7);

index = 1;
for clusterIndex = 1:numberOfClusters
    
   %%% number of subpaths
   numberOfSubpathComponents = numberOfClusterSubPaths(clusterIndex);
   
   %%% subpath delays in the same time cluster
   subpathDelays = t_mn.(['c',num2str(clusterIndex)]);
   subpathPowers = subPathPowers.(['c',num2str(clusterIndex)]);
   subpathPhases = phases_mn.(['c',num2str(clusterIndex)]); 
   
   for subpathIndex = 1:numberOfSubpathComponents
       
       %%% extract subpath delay, power, phase
       subpathDelay = subpathDelays(subpathIndex);
       subpathPower = subpathPowers(subpathIndex)*10^(Pr_dBm/10);
       subpathPhase = subpathPhases(subpathIndex);
       
       %%% extract AOD Azi/El
       subpath_AOD_AziEL = subpath_AODs.(['c',num2str(clusterIndex)]).(['SP',num2str(subpathIndex)]);

       %%% extract AOA Azi/El
       subpath_AOA_AziEL = subpath_AOAs.(['c',num2str(clusterIndex)]).(['SP',num2str(subpathIndex)]);
   
       %%% rearrange in row vector
       storeMat = [subpathDelay subpathPower subpathPhase subpath_AOD_AziEL subpath_AOA_AziEL];
       
       %%% store
       powerSpectrum(index,:) = storeMat;
       
       %%% increment for next cycle
       index = index+1;       
       
   end%%end of subpathIndex
    
end%%end of clusterIndex


end