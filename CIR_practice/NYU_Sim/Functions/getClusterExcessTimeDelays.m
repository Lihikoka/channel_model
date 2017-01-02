function tau_n = getClusterExcessTimeDelays(mu_tau,rho_mn,minVoidInterval)
% Generate time cluster delays.
%
% Inputs:
%   - mu_tau: mean excess delay in ns
%   - rho_mn: structure of intra-cluster subpath delays in ns
%   - minVoidInterval: minimum inter-cluster void interval, typically set
%   to 25 ns for outdoor environments
% Output:
%   - tau_n: array of time cluster excess delays in ns
%
% Copyright © 2016 NYU

    %%% number of time clusters
    numberOfTimeClusters = size(fieldnames(rho_mn),1);

    %%% initialize cluster delays
    tau_n = zeros(1,numberOfTimeClusters);
        
    %%% void interval between consecutive clusters
    clusterVoidInterval = minVoidInterval+2.5;
    
    %%%% generate cluster delays as exponential
    tau_n_prime = exprnd(mu_tau,[1 numberOfTimeClusters]);
        
    %%% normalize
    tau_n_double_prime = sort(tau_n_prime - min(tau_n_prime));
    
    temp = rho_mn.c1(end);
    %%% this for loop starts at the 2nd cluster index because the first
    %%% cluster index is always 0
    for clusterIndex = 2:numberOfTimeClusters        
        
        %%% add the last cluster sub-path delay of the previous cluster to
        %%% the current cluster delay for no overlap in multipath
        %%% components. 
        
        tau_n(clusterIndex) = tau_n_double_prime(clusterIndex)+temp+clusterVoidInterval;        
        
        %%% cluster sub-path delays of the previous cluster
        rho_m = rho_mn.(['c',num2str(clusterIndex)]);
        
        %%% keep track of the last intra-cluster delay
        temp = tau_n(clusterIndex)+rho_m(end);
        
    end
    
        

end




