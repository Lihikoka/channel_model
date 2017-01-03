clc
clear

Nt = 81; Nr = 1;  % Number of antenna at Tx and Rx
K = 20; % Number of user
L = 2; % Number of path

phi = rand(K, L)*pi - pi/2; % K users' azimuth in radian in [-pi/2, pi/2]
phi*180/pi
theta = 0.5*sin(phi);
U = zeros(Nt, Nt);
I = zeros(Nt, 1);
a = zeros(Nt, L);

% Parameters for path loss model
% 28GHz, TX height=7, Rx height=1.5
alpha = 75.85; % in dB
beta_bar = 3.73;
d = 100; % TRSep
PL_dB = alpha + 3.73*10*log(d);
PL = 10^(PL_dB/10);

% Beam index set
for l=1:Nt
    I(l) = (l-1) - (Nt-1)/2;
end

% Channel matrix
H_LOS = [];
H = [];
for k = 1:K
    % beta = (randn(1)+1i*randn(1))/sqrt(2); % Rayleigh fading complex path loss
    beta(1) = exp(-1i*2*pi*rand(1)); % |beta|^2=1
    beta(2) = 0.1*exp(-1i*2*pi*rand(1));
    for n = 1:Nt
        for l = 1:L
            a(n, l) = exp(-1i*2*pi*theta(k,l)*I(n));
        end
    end
    h_k_LOS = beta(1)*a(:, 1);
    H_LOS = horzcat(H_LOS, h_k_LOS);
    h_k = h_k_LOS;
    for i = 2:L
        h_k = h_k + beta(i)*a(:, i);
    end
    H = horzcat(H, h_k);
end

% DFT matrix
for i = 1:Nt
    for j = 1:Nt
        U(i,j) = exp(-1i*2*pi*I(j)/Nt*I(i))/sqrt(Nt);
    end
end

% Beamspace channel
H_b_LOS = ctranspose(U)*H_LOS; % Only LOS
H_b_LOS_power = abs(ctranspose(H_b_LOS)).^2;

H_b = ctranspose(U)*H; % LOS and 1 NLOS
H_b_power = abs(ctranspose(H_b)).^2;

figure(1)
contour(-90/Nt+180/Nt*I, 1:K, H_b_LOS_power);
xlabel('TX BEAM DIRECTION');
ylabel('MOBILE STATION');

figure(2)
contour(-90/Nt+180/Nt*I, 1:K, H_b_power);
xlabel('TX BEAM DIRECTION');
ylabel('MOBILE STATION');

% for i = 1:K
%     polarplot(-pi/2:pi/(Nt-1):pi/2, H_b_LOS_power(i, :));
%     hold on
% %     stem(I, H_b_power(i, :));
% end
%rlim([0.5 0])





%%%%%%%%%%%%% NYU_Sim %%%%%%%%%%%%% 
% frequency = '2873_GHz';
% scenario = 'outdoor';
% N=40; % N users
% K_dB = 5; % 5~15
% 
% string = ['load Results_',frequency,'_',scenario,'_',num2str(N),'_CIRs'];
% eval(string);
% 
% 
% indexLOS = [];
% % fetch CIR data
% for i = 1:N
%     string = ['CIR(', num2str(i), ') = ', 'CIR_Struct.CIR_', num2str(i), ';'];
%     eval(string);
%     % allPower = [allPower; CIR(i).pathPowers];
%     stem(CIR(i).pathDelays, CIR(i).pathPowers);
%     if CIR(i).pathPowers(1) >= max(CIR(i).pathPowers)
%         indexLOS = [indexLOS, i];
%     end
%     numTx = CIR(i).NumOfTxElements;
%     numRx = CIR(i).NumOfRxElements;
% end
% 
% % Simulation
% % N = length(indexLOS);
% % Tx_power = 1;
% % noise_power = 0.001; % sigma^2
% % H = zeros(numTx, N);
% % H = CIR(1).H{1};
% H = [];
% for i=1:N
%     H = vertcat(H, CIR(i).H{1,1}); % Users' first path channel matrix
% end
% % DLA
% U = zeros(numTx, numTx);
% I = zeros(numTx, 1);
% for l=1:numTx
%     I(l) = (l-1) - (numTx-1)/2; % Beam index set
% end
% for i=1:numTx
%     for j=1:numTx
%         U(i,j) = exp(-1i*2*pi*I(j)/numTx*I(i))/sqrt(numTx); % DFT matrix
%     end
% end
% 
% H_b = H*U; % Beamspace channel matrix
% s = (randn(N*numRx,1) + 1i*randn(N*numRx,1))/sqrt(2); % symbol vertor
% 
% % Compute beamspace channel power
% for i=1:size(H_b, 2)
%     H_power = conj(H).*H;
%     H_b_power =conj(H_b).*H_b;
% end
% 
% plot(I, transpose(conj(H_b).*H_b));
% 
% % Normalization
% for i=1:size(H_b, 1)
%     H_b_power_max = 0;
%     H_power_max = 0;
%     for j=1:size(H_b, 2)
%         if H_b_power(i,j) > H_b_power_max
%             H_b_power_max = H_b_power(i,j);
%         end
%         if H_power(i,j) > H_power_max
%             H_power_max = H_power(i,j);
%         end
%     end
%     H_power(i,:) = H_power(i,:)/H_power_max;
%     H_b_power(i,:) = H_b_power(i,:)/H_b_power_max;
% end
% 
% % contour(H_b_power);
% % stem(I, H_b_power);
% for i=1:length(indexLOS)
% % for i=1:N
% %     figure(1)
% %     stem(I, H_power(indexLOS(i),:))
% %     hold on
%     figure(i)
%     stem(I, H_b_power(indexLOS(i),:))
% %     stem(I, H_b_power(i,:))
%     hold off
%     legend('beamspace channel power');
% end
% 
% xlabel('TX BEAM index');
% ylabel('|H_b^H|^2');


% Precoding Matrix
% for i=1:N
%     F_MF = H; % Matched filter
%     F_ZF = H/(ctranspose(H)*H); % Zero-forcing 
%     eta = noise_power*numRx/Tx_power; % Regulization factor pf wiener filter
%     F_WF = (H*ctranspose(H)+eta*eye(numRx))\H; % Wiener filter
% 
%     alpha_MF = sqrt(Tx_power/trace(F_MF*ctranspose(F_MF)));
%     alpha_ZF = sqrt(Tx_power/trace(F_ZF*ctranspose(F_ZF)));
%     alpha_WF = sqrt(Tx_power/trace(F_WF*ctranspose(F_WF)));
% 
%     G_MF = alpha_MF * F_MF;
%     G_ZF = alpha_ZF * F_ZF;
%     G_WF = alpha_WF * F_WF;
%     
%     % Beamspace channel matrix and precoding matrix
%     
%     G_b_MF = U\G_MF;
%     G_b_ZF = U\G_ZF;
%     G_b_WF = U\G_WF;
%     
%     r_MF = ctranspose(H)*G_MF*s;
%     r_b_MF = ctranspose(H_b)*(G_b_MF)*s; %+ noise_power*(randn(N,1)+1i*randn(N,1))/sqrt(2);
% end







%%%%%%%%% generate Rician Fading channel matrix

% Generate Uncorrelated Rician Fading channel matrix
% for i=1:1
%     L = length(CIR(indexLOS(i)).pathDelays);
%     A = Ric_channel_matrix(numTx, numRx, K_dB, L);
% end
% 
% % Generate Tx and Rx correlation matrix
% Rtx = zeros(numTx, numTx);
% Rrx = zeros(numRx, numRx);
% for i=1:numRx
%     for j=i:numRx
%         %Rrx(i, j) = 
%     end
% end
% 
% 
% 

% %%%%%%%%%% plot AoA distribution
% AOARange = 0:0.01:360;
% AOAGraph = zeros(1, 36001);
% AOApdf = zeros(1, 36001);
% 
% for i = 1:length(allAOA)
%     pos          = find(AOARange>=allAOA(i)); 
% %     AOAIndex = round(AOARange(pos(1)));
%     AOAIndex = ((AOARange(pos(1))-allAOA(i))*(pos(1)-1) + ... % interpolation
%                 (allAOA(i)-AOARange(pos(1)-1))*(pos(1))) / ...
%                 (AOARange(pos(1))-AOARange(pos(1)-1));
%     AOAIndex = round(AOAIndex);
%     AOAGraph(AOAIndex) = AOAGraph(AOAIndex) + 1;
% end
% 
% sum = 0;
% for i = 1:length(AOARange)
%     sum = sum + AOAGraph(i);
%     if ~mod(AOARange(i), 1)
%         AOApdf(i) = sum/length(allAOA);
%         sum=0;
%     end
% end
% 
% 
% figure(1)
% plot(AOARange, AOApdf)
% xlabel('AOA(deg)');
% ylabel('number');
% 
