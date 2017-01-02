% Copyright © 2016 NYU

% Reference distance d_0 (m)
d0 = 1;

% Applicable distance ranges for the channel model
dmin = 30; dmax = 60;
        
% path loss parameters - SF in dB
n = 2; SF = 3.6;

% mean number of AOD spatial lobes
mu_AOD = 1.9;

% mean number of AOA spatial lobes
mu_AOA = 1.8;

% X_max
X_max = 0.2;

% Mean excess delay mu_tau (ns)
mu_tau = 123;

% Minimum inter-cluster void interval (ns)
minVoidInterval = 25;

% Per-cluster shadowing (dB)
sigmaCluster = 1;

% Time cluster decay constant (ns)
Gamma = 25.9;

% Per-subpath shadowing (dB)
sigmaSubpath = 6;

% Subpath decay constant (ns)
gamma = 16.9;

% Mean zenith of departure (ZOD) (degrees)
mean_ZOD = -12.6;

% Standard deviation of ZODs (degrees)
sigma_ZOD = 5.9;

% Azimuth offset parameter (degrees) from AOD spatial lobe centroid
std_AOD_RMSLobeAzimuthSpread = 8.5;

% Elevation offset parameter (degrees) from ZOD spatial lobe centroid
std_AOD_RMSLobeElevationSpread = 2.5;

% Distribution type that generates the elevation offsets from AOD spatial
% lobe centroid
distributionType_AOD = 'Gaussian'; 

% Mean zenith of arrival (ZOA) (degrees)
mean_ZOA = 10.8;

% Standard deviation of ZOAs (degrees)
sigma_ZOA = 5.3;

% Azimuth offset parameter (degrees) from AOA spatial lobe centroid
std_AOA_RMSLobeAzimuthSpread = 10.5;

% Elevation offset parameter (degrees) from ZOA spatial lobe centroid
std_AOA_RMSLobeElevationSpread = 11.5;

% Distribution type that generates the elevation offsets from AOA spatial
% lobe centroid
distributionType_AOA = 'Laplacian'; 













