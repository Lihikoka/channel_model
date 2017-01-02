% Copyright © 2016 NYU

% Reference distance d_0 (m)
d0 = 1;

% Applicable distance ranges for the channel model
dmin = 60; dmax = 200;
        
% path loss parameters - SF in dB
n = 3.3; SF = 7.6;

% mean number of AOD spatial lobes
mu_AOD = 1.5;

% mean number of AOA spatial lobes
mu_AOA = 2.5;

% X_max
X_max = 0.5;

% Mean excess delay mu_tau (ns)
mu_tau = 83;

% Minimum inter-cluster void interval (ns)
minVoidInterval = 25;

% Per-cluster shadowing (dB)
sigmaCluster = 3;

% Time cluster decay constant (ns)
Gamma = 56;

% Per-subpath shadowing (dB)
sigmaSubpath = 6;

% Subpath decay constant (ns)
gamma = 15.3;

% Mean zenith of departure (ZOD) (degrees)
mean_ZOD = -4.9;

% Standard deviation of ZODs (degrees)
sigma_ZOD = 4.5;

% Azimuth offset parameter (degrees) from AOD spatial lobe centroid
std_AOD_RMSLobeAzimuthSpread = 7;

% Elevation offset parameter (degrees) from ZOD spatial lobe centroid
std_AOD_RMSLobeElevationSpread = 3.5;

% Distribution type that generates the elevation offsets from AOD spatial
% lobe centroid
distributionType_AOD = 'Gaussian'; 

% Mean zenith of arrival (ZOA) (degrees)
mean_ZOA = 3.6;

% Standard deviation of ZOAs (degrees)
sigma_ZOA = 4.8;

% Azimuth offset parameter (degrees) from AOA spatial lobe centroid
std_AOA_RMSLobeAzimuthSpread = 6;

% Elevation offset parameter (degrees) from ZOA spatial lobe centroid
std_AOA_RMSLobeElevationSpread = 3.5;

% Distribution type that generates the elevation offsets from AOA spatial
% lobe centroid
distributionType_AOA = 'Laplacian'; 



