function numberOfClusterSubPaths = getNumberOfClusterSubPaths(numberOfTimeClusters)
% Generate the number of intra-cluster subpaths
%
% Inputs:
%   - numberOfTimeClusters: the number of time clusters
% Output:
%   - numberOfClusterSubPaths: an array containing the number of subpaths
%   for each time cluster
%
% Copyright © 2016 NYU

    % Generate number of subpaths in each time cluster
    numberOfClusterSubPaths = randi([1 30],[numberOfTimeClusters 1]);


end