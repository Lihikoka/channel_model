function rho_mn = getIntraClusterDelays(numberOfClusterSubPaths,X_max)
% Generate intra-cluster subpath delays in ns.
%
% Inputs:
%   - numberOfClusterSubPaths: array containing the number of subpaths for
%   each time cluster
%   - X_max: a number between 0 and 1
% Output:
%   - rho_mn: a structure containing intra-cluster subpath delays, in ns
%
% Copyright © 2016 NYU



    %%% initialize the structure that will contain the intra cluster delays
    rho_mn = struct;
    
    %%% number of clusters
    numberOfClusters = size(numberOfClusterSubPaths,1);

    %%% for loop iterates N times for each cluster
    for clusterIndex=1:numberOfClusters

        %%% number of sub-paths in current cluster
        numberOfComponents = numberOfClusterSubPaths(clusterIndex);

        %%% generate a set of component delays
        arrayTemp = 2.5*(1:numberOfComponents);

        %%% field name
        str = ['c',num2str(clusterIndex)];

        %%% sort the array
        sortedArray = sort(arrayTemp - min(arrayTemp));
        
        %%% store the components
        X = X_max*rand;
        rho_mn.(str) = sortedArray.^(1+X);        

    end

end