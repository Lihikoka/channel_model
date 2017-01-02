function [numberOfTimeClusters,numberOfAOALobes,numberOfAODLobes] = ...
    getNumClusters_AOA_AOD(mu_AOA,mu_AOD)  
% Generate the number of time clusters, the number of AOD spatial lobes,
% and the number of AOA spatial lobes 
%
% Inputs:
%   - mu_AOA: the mean number of spatial lobes at the RX 
%   - mu_AOD: the mean number of spatial lobes at the TX 
% Output:
%   - numberOfTimeClusters: the number of time clusters
%   - numberOfAOALobes: the number of spatial lobes at the RX
%   - numberOfAODLobes: the number of spatial lobes at the TX
%
% Copyright © 2016 NYU

     numberOfTimeClusters = randi([1,6]);
     
     aod_instance = poissrnd(mu_AOD);
     
     numberOfAODLobes = max(1,min(5,aod_instance));
     
     aoa_instance = poissrnd(mu_AOA);
     
     numberOfAOALobes = max(1,min(5,aoa_instance));
     
end