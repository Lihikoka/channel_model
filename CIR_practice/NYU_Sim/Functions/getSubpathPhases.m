function phases_mn = getSubpathPhases(rho_mn)

%%% initialize struct
phases_mn = struct;

%%% number of time clusters
numberOfTimeClusters = size(fieldnames(rho_mn),1);

for clusterIndex = 1:numberOfTimeClusters
    
    %%% the intra-cluster subpath delays
    rho_m = rho_mn.(['c',num2str(clusterIndex)]);
    
    %%% number of subpaths
    numberOfSubpaths = numel(rho_m);
    
    %%% initialize subpath phases
    subpathPhases = 2*pi*rand(1,numberOfSubpaths);
        
    %%% store subpath phases
    phases_mn.(['c',num2str(clusterIndex)]) = subpathPhases;
    
end%%end of clusterIndex for loop


end