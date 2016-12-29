% NYU-Simulator Version 1.1 (NYU-Sim), developed by:
%
% Mathew K. Samimi, Shu Sun - NYU WIRELESS, 16 February 2016
%
% This script runs the NYU-Simulator
%
% Variable inputs to this script:
%   N (number of CIRs to generate)
%   See loadUserInputs.m in the ReadFiles folder
%
% Outputs of this script:
% ['Results_',frequency,'_',environment,'_',scenario,'_',num2str(N),'_CIRs.mat']
% ['SecondOrderStats_',frequency,'_',environment,'_',scenario,'_',num2str(N),'_CIRs.txt'];
%
% Copyright ?2016 NYU
%
% See REFERENCES below for more information:
%
% [1] M. Samimi et al., "28 GHz Angle of Arrival and Angle of Departure 
% Analysis for Outdoor Cellular Communications Using Steerable Beam 
% Antennas in New York City," in 2013 IEEE  Vehicular Technology 
% Conference (VTC Spring), pp.1-6, 2-5 June 2013.
% URL: http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=6691812&isnumber=6691801
% 
% [2] M. K. Samimi, T. S. Rappaport, "Ultra-wideband statistical channel
% model for non line of sight millimeter-wave urban channels," in 2014
% IEEE Global Communications Conference (GLOBECOM), pp. 3483-3489, 8-12
% Dec. 2014.
% URL: http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=7037347&isnumber=7036769
% 
% [3] M. K. Samimi, T. S. Rappaport, "3-D statistical channel model for
% millimeter-wave outdoor mobile broadband communications," 2015 IEEE 
% International Conference on in Communications (ICC), pp.2430-2436, 
% 8-12 June 2015.
% URL: http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=7248689&isnumber=7248285
% 
% [4] M. K. Samimi, T. S. Rappaport, “Statistical Channel Model with 
% Multi-Frequency and Arbitrary Antenna Beamwidth for Millimeter-Wave
% Outdoor Communications,?in 2015 IEEE Global Communications Conference,
% Exhibition & Industry Forum (GLOBECOM) Workshop, Dec. 6-10, 2015.
% URL: http://arxiv.org/abs/1510.03081
% 
% [5] M. K. Samimi, T. S. Rappaport, ?8 GHz Millimeter-Wave Ultrawideband 
% Small-Scale Fading Models in Wireless Channels,?in 2016 IEEE 
% Vehicular Technology Conference (VTC2016-Spring), 15-18 May, 2016.
% URL: http://arxiv.org/abs/1511.06938
% 
% [6] M. K. Samimi, S. Sun, and T. S. Rappaport, “MIMO Channel Modeling 
% and Capacity Analysis for 5G Millimeter-Wave Wireless Systems,?in 
% the 10th European Conference on Antennas and Propagation (EuCAP?016), 
% April 2016.
% URL: http://arxiv.org/abs/1511.06940
% 
% [7] M. K. Samimi, T. S. Rappaport, “Local Multipath Model Parameters for 
% Generating 5G Millimeter-Wave 3GPP-like Channel Impulse Response,?
% in the 10th European Conference on Antennas and Propagation (EuCAP?016),
% April 2016.
% URL: http://arxiv.org/abs/1511.06941
 
clear 
clc
close all

% master folder, and add path to all underlying subfolders
%setMasterFolder

% Folder where the output is saved;
resultsFolder = ['D:\Users\yenhua\Documents\github\channel_model\CIR_practice\'];

% Load user inputs
%loadUserInputs

%%%%%%%%%%%%%%%%%
% NYU-Simulator %
%%%%%%%%%%%%%%%%%

% Generate N Channel Impulse Reponses (CIRs): one for each of the N users
% Output of generateCIR: CIR_Struct
% generateCIR

% Output the generated CIRs
% CIR_Struct;

% Generate the local area CIRs based on the previously generated single
% CIRs
% CIR_Struct = getLocalCIR(CIR_Struct,TxArrayType,RxArrayType,Nt,Nr,Wt,Wr,dTxAnt,dRxAnt);

% Compute RMS DS and RMS Angular Spreads, and output to .txt files
% cd(resultsFolder)
%  AS = extractStatistics(CIR_Struct)

% Save output for later processing
% save(['Results_',frequency,'_',scenario,'_',num2str(N),'_CIRs'],'CIR_Struct')
frequency = '28_GHz';
scenario = 'outdoor';
N=100;
load([resultsFolder, 'Results_', frequency,'_',scenario,'_',num2str(N),'_CIRs'])

%===
% Plotting
% h1_delay = CIR_Struct.CIR_1.pathDelays;
% h1_power = CIR_Struct.CIR_1.pathPowers;
% 
% h2_delay = CIR_Struct.CIR_2.pathDelays;
% h2_power = CIR_Struct.CIR_2.pathPowers;
% 
% subplot(1,2,1);
% stem(h1_delay, h1_power, 'r');
% xlabel('Path delay (ns)');
% ylabel('Path power (mW)');
% 
% subplot(1,2,2);
% stem(h2_delay, h2_power, 'b');
% xlabel('Path delay (ns)');
% ylabel('Path power (mW)');


H_delay=[];
H_aoa=[];
H_aod=[];
H_power=[];
for ii = 1:N
    H_delay= cat(2,H_delay, CIR_Struct.(['CIR_', num2str(ii)]).pathDelays');
    H_aoa= cat(2,H_aoa, CIR_Struct.(['CIR_', num2str(ii)]).AOAs');
    H_aod= cat(2,H_aod, CIR_Struct.(['CIR_', num2str(ii)]).AODs');
    H_power= cat(2,H_power, CIR_Struct.(['CIR_', num2str(ii)]).pathPowers');
end
% set(gca,'xTick',-2*pi:pi/2:2*pi);
% set(gca,'xTickLabel',{'0','2*pi','-2*pi'});
H_aoa_wrap = wrap(H_aoa,1);
% [Asorted,I_a]=sort(H_aoa);
% plot(Asorted,(I_a));
% H_aoa = mod(H_aoa,360);

 fh1 = figure(1);
% [data,bin] = ksdensity(H_aoa_wrap,'bandwidth',0.2,'function','pdf','kernel','epanechnikov');
[data,bin] = ksdensity(H_power,'bandwidth',0.1,'function','pdf','kernel','epanechnikov','npoints',length(H_power));
plot(bin,data,'r-', 'LineWidth',1.5);
% nbins = 25;
% histogram(H_aoa,nbins,'Normalization','pdf');

% xx = hist(H_aoa_wrap,[-2*pi:2*pi]);
% plot([-2*pi,2*pi],xx./numel(xx));

hold on;
% sigma = std(H_aoa_wrap); %5.3*pi/180;
% sigma = std(H_power);
% sigma = AS;
sigma = 4.8*pi/180;

%Truncated Laplacian
t = -pi:0.1:pi;
% t = -360:1:360;
PAS = zeros(1,length(t));
Q= 1/( 1-exp(-sqrt(2)*pi/sigma) );
phi_0 = 0*pi/180; % mean
for ii=1:length(t)
    if t(ii) > phi_0-pi & t(ii) < phi_0+pi
%     if t(ii) > 180 & t(ii) < 180
       PAS(ii) = Q/sqrt(2)/sigma*exp( -sqrt(2)*abs(t(ii)-phi_0)/sigma );
    end
end
plot(t,PAS,'k-.','LineWidth',1.5);

% Uniform PAS
% d = -2*pi:0.1:2*pi;
% DeltaPhi = sqrt(3)*sigma;
% Q = 1/2/DeltaPhi;
% PAS = zeros(1,length(t));
% for ii=1:length(t)
%     if t(ii) > -DeltaPhi & t(ii) < DeltaPhi
%        PAS(ii) = Q;
%     end
% end
% plot(t,PAS,'b--','LineWidth', 1.5);

xlhand = get(gca,'xlabel');
set(xlhand,'string','Angle (radian)','fontsize',14);

ylhand = get(gca,'ylabel');
set(ylhand,'string','PAS','fontsize',14);

legend('Simulated','Laplace');
% legend('Simulated','Laplace','Uniform');
% filename = ['PAS_frequency',frequency,'_',scenario,'_',num2str(N),'_CIRs'];
% saveas(fh1, filename, 'epsc');

% fh2 = figure(2);
% set(fh2, 'color', 'white');
% xi=linspace(min(H_delay),max(H_delay),100);
% yi=linspace(min(H_aoa),max(H_aoa),100);
% [XI YI]=meshgrid(xi,yi);
% ZI = griddata(H_delay,H_aoa,H_power,XI,YI);
% mesh(XI,YI,ZI);
% 
% xlhand = get(gca,'xlabel');
% set(xlhand,'string','Delay (ns)','fontsize',14);
% % xlim([100 1000]);
% 
% ylhand = get(gca,'ylabel');
% set(ylhand,'string','AoA (degree)','fontsize',14);
% 
% zlhand = get(gca,'zlabel');
% set(zlhand,'string','Power (mW)','fontsize',14);
% filename = ['Delay-AOA_',frequency,'_',scenario,'_',num2str(N),'_CIRs'];
% saveas(fh2, filename, 'epsc');

% fh3 = figure(3);
% set(fh3, 'color', 'white');
% xi=linspace(min(H_delay),max(H_delay),100);
% yi=linspace(min(H_aod),max(H_aod),100);
% [XI YI]=meshgrid(xi,yi);
% ZI = griddata(H_delay,H_aod,H_power,XI,YI);
% mesh(XI,YI,ZI);
% 
% xlhand = get(gca,'xlabel');
% set(xlhand,'string','Delay (ns)','fontsize',14);
% % xlim([100 1000]);
% 
% ylhand = get(gca,'ylabel');
% set(ylhand,'string','AoD (degree)','fontsize',14);
% 
% zlhand = get(gca,'zlabel');
% set(zlhand,'string','Power (mW)','fontsize',14);
% filename = ['Delay-AOD_',frequency,'_',scenario,'_',num2str(N),'_CIRs'];
% saveas(fh3, filename, 'epsc');














