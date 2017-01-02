function [angularInfo_struct, cluster_subpath_lobe_mapping] = ...
    getSubpathAngles(numberOfLobes,numberOfClusterSubPaths,mean_zenith,...
    sigma_zenith,std_RMSLobeElevationSpread,std_RMSLobeAzimuthSpread,...
    distributionType)

% Generate the multipath azimuth and elevation angles
%
% Inputs:
%   - numberOfLobes: the number of spatial lobes
%   - numberOfClusterSubPaths: an array containing the number of subpaths
%   in each time cluster
%   - mean_zenith: mean zenith angle, in degrees
%   - sigma_zenith: standard deviation of the zenith angle distribution, in
%   degrees
%   - std_RMSLobeElevationSpread: standard deviation of the elevation
%   offset from the lobe centroid
%   - std_RMSLobeAzimuthSpread: standard deviation of the azimuth offset
%   from the lobe centroid
%   - distributionType: a string specifying which distribution to use:
%   'Gaussian' or 'Laplacian'
% Output:
%   - angularInfo_struct: a structure containing the subpath angles
%   - cluster_subpath_lobe_mapping: an array containing a mapping between
%   cluster number, subpath number, and lobe number, for each multipath
%
% Copyright © 2016 NYU

%%% initialize struct
angularInfo_struct = struct;

%%% generate mean lobe angles
%%% lobe indices
lobeIndices = 1:numberOfLobes;

%%% discretize the azimuth plane in the same number of lobes
theta_min_array = 360*(lobeIndices-1)/numberOfLobes;
theta_max_array =  360*lobeIndices/numberOfLobes;

%%% initialize the lobe aoas to 0
mean_AzimuthAngles = zeros(numberOfLobes,1);
mean_ElevationAngles = mean_zenith+sigma_zenith*randn(numberOfLobes,1);

%%% generate a lobe aoa that is uniformly random in each section of the
%%% azimuth plane
for lobeIndex=1:numberOfLobes
    %%% azimuth angle
    aziAngle = theta_min_array(lobeIndex)+(theta_max_array(lobeIndex)-theta_min_array(lobeIndex))*rand;
    mean_AzimuthAngles(lobeIndex) = aziAngle;
end%%end of lobeIndex

numberOfClusters = numel(numberOfClusterSubPaths);

%%% total number of subpath components
totalNumberOfSubpathComponents = sum(numberOfClusterSubPaths);

%%% initialize the mapping matrix
cluster_subpath_lobe_mapping = zeros(totalNumberOfSubpathComponents,3);

index = 1;
%%% generate subpath angles
for clusterIndex = 1:numberOfClusters    
       
    %%% number of components
    numberOfClusterComponents = numberOfClusterSubPaths(clusterIndex);
    
    for subpathIndex = 1:numberOfClusterComponents
        
        %%% assign random lobe index
        randomLobeIndex = randi([1 numberOfLobes]);
        
        %%% extract mean azi/el angles
        meanAzi = mean_AzimuthAngles(randomLobeIndex);
        meanEl = mean_ElevationAngles(randomLobeIndex);

        %%% Generate azimuth offsets from spatial lobe centroid
        deltaAzi = std_RMSLobeAzimuthSpread*randn;
        
        %%% generate elevation offsets within spatial lobes
        switch distributionType
            case 'Gaussian'
                deltaEl = std_RMSLobeElevationSpread*randn;                
            case 'Laplacian'
                z = -.5+rand;
                b = std_RMSLobeElevationSpread/sqrt(2);
                deltaEl = -b*sign(z).*log(1-2*abs(z));                
            otherwise
        end
        
        %%% compute subpath angles
        subpathAzi = mod(meanAzi+deltaAzi,360);
        subpathEl = min(max(meanEl+deltaEl,-60),60);                        
        
        %%% store
        angularInfo_struct.(['c',num2str(clusterIndex)]).(['SP',num2str(subpathIndex)])...
            = [subpathAzi subpathEl];        
        
        %%% store mapping
        cluster_subpath_lobe_mapping(index,:) = [clusterIndex subpathIndex randomLobeIndex];
        index = index+1;
        
    end%%end of subpathIndex
    
end%%end of clusterIndex



end