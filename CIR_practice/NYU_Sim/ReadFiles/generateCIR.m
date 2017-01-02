% Copyright ?2016 NYU

% Output to command window
disp(['Generating ',num2str(N),' CIRs...'])

% Allocate number of LOS and NLOS users
N_LOS = round(N*PercentLOS);
N_NLOS = N-N_LOS;

% structure containing generated CIRs
CIR_Struct = struct;

% Generate the users, according to LOS or NLOS
CIRStoreIdx = 1;
for CIRIdx = 1:N_LOS    
   CIR_Struct.(['CIR_',num2str(CIRStoreIdx)]).environment = 'LOS';
   CIRStoreIdx = CIRStoreIdx+1;
end
for CIRIdx = 1:N_NLOS
   CIR_Struct.(['CIR_',num2str(CIRStoreIdx)]).environment = 'NLOS'; 
   CIRStoreIdx = CIRStoreIdx+1;
end

% Generate the N CIRs
for CIRIdx = 1:N
    
    % Extract environment (LOS or NLOS)
    environment = CIR_Struct.(['CIR_',num2str(CIRIdx)]).environment; 
    
    % Load desired parameters
    eval([scenario,'_',environment,'_',accessType,'_',frequency,'_ChannelParams'])
    
    % Step 1: Generate T-R Separation distance (m) ranging from dmin - dmax.
    TRDistance = getTRSep(dmin,dmax);
    
    % Step 2: Generate the total received omnidirectional power (dBm) and 
    % path loss (dB) 
    [Pr_dBm, PL_dB]= getRXPower(frequency,n,SF,TXPower,TRDistance,d0,dynamicRange);

    % Step 3: Generate the number of time clusters N, and number of AOD and
    % AOA spatial lobes
    [numberOfTimeClusters,numberOfAOALobes,numberOfAODLobes] = ...
        getNumClusters_AOA_AOD(mu_AOA,mu_AOD);

    % Step 4: Generate the number of cluster subpaths M_n for each time
    % cluster
    numberOfClusterSubPaths = getNumberOfClusterSubPaths(numberOfTimeClusters);

    % Step 5: Generate the intra-cluster subpath delays rho_mn (ns)
    rho_mn = getIntraClusterDelays(numberOfClusterSubPaths,X_max);

    % Step 6: Generate the phases (rad) for each cluster
    phases_mn = getSubpathPhases(rho_mn);

    % Step 7: Generate the cluster excess time delays tau_n (ns)
    tau_n = getClusterExcessTimeDelays(mu_tau,rho_mn,minVoidInterval);

    % Step 8: Generate temporal cluster powers (normalized to 1 mW)
    clusterPowers = getClusterPowers(tau_n,Gamma,sigmaCluster);

    % Step 9: Generate the cluster subpath powers (normalized to 1 mW)
    subpathPowers = ...
        getSubpathPowers(rho_mn,clusterPowers,gamma,sigmaSubpath,TXPower,dynamicRange);

    % Step 10: Recover absolute propagation times t_mn (ns) of each subpath 
    % component
    t_mn = getAbsolutePropTimes(TRDistance,tau_n,rho_mn);

    % Step 11: Recover AODs and AOAs of the multipath components
    [subpath_AODs, cluster_subpath_AODlobe_mapping] = ...
        getSubpathAngles(numberOfAODLobes,numberOfClusterSubPaths,mean_ZOD,...
        sigma_ZOD,std_AOD_RMSLobeElevationSpread,std_AOD_RMSLobeAzimuthSpread,...
        distributionType_AOD);
    [subpath_AOAs, cluster_subpath_AOAlobe_mapping] = ...
        getSubpathAngles(numberOfAOALobes,numberOfClusterSubPaths,mean_ZOA,...
        sigma_ZOA,std_AOA_RMSLobeElevationSpread,std_AOA_RMSLobeAzimuthSpread,...
        distributionType_AOA);

    % Step 12: Construct the multipath parameters
    powerSpectrum = getPowerSpectrum(numberOfClusterSubPaths,t_mn,subpathPowers,phases_mn,...
        subpath_AODs,subpath_AOAs,Pr_dBm);

    % Construct the 3-D lobe power spectra at TX and RX
    AOD_LobePowerSpectrum = getLobePowerSpectrum(numberOfAODLobes,cluster_subpath_AODlobe_mapping,powerSpectrum,'AOD');
    AOA_LobePowerSpectrum = getLobePowerSpectrum(numberOfAOALobes,cluster_subpath_AOAlobe_mapping,powerSpectrum,'AOA');

    % Store CIR parameters
    CIR.pathDelays = powerSpectrum(:,1);
    CIR.pathPowers = powerSpectrum(:,2);
    CIR.pathPhases = powerSpectrum(:,3);
    CIR.AODs = powerSpectrum(:,4);
    CIR.ZODs = powerSpectrum(:,5);
    CIR.AOAs = powerSpectrum(:,6);
    CIR.ZOAs = powerSpectrum(:,7);

    % Various other information for this CIR
    CIR.frequency = frequency;
    CIR.TXPower = TXPower;
    CIR.OmniPower = Pr_dBm;
    CIR.OmniPL = PL_dB;
    CIR.TRSep = TRDistance;
    CIR.environment = environment;
    CIR.scenario = scenario;
    CIR.accessType = accessType;
    CIR.HPBW_TX = [theta_3dB_TX phi_3dB_TX];
    CIR.HPBW_RX = [theta_3dB_RX phi_3dB_RX];
    
    % Store
    CIR_Struct.(['CIR_',num2str(CIRIdx)]) = CIR;

end% end of CIRIdx

% Output to command window
disp('Generating CIRs..Done.')

