clc
clear

frequency = '2873_GHz';
scenario = 'outdoor';
N=20; % N users
K_dB = 5; % 5~15

string = ['load Results_',frequency,'_',scenario,'_',num2str(N),'_CIRs'];
eval(string);

allAOA = [];
allPower = [];
indexLOS = [];

% fetch CIR data
for i = 1:N
    string = ['CIR(', num2str(i), ') = ', 'CIR_Struct.CIR_', num2str(i), ';'];
    eval(string);
    allAOA = [allAOA; CIR(i).AOAs];
    allPower = [allPower; CIR(i).pathPowers];
    if strcmp(CIR(i).environment, 'LOS')
        indexLOS = [indexLOS, i];
    end
    numTx = CIR(i).NumOfTxElements;
    numRx = CIR(i).NumOfRxElements;
    for j=1:length(CIR(i).H)
        H = CIR(i).H(j);
    end
end

% Simulation
% N = length(indexLOS);
Tx_power = 1;
noise_power = 0.001; % sigma^2
H = zeros(numTx, N);
for i=1:N
    H(:,i) = CIR(i).H{1}; % 25 LOS users' first path CIR
end
% DLA: DFT matrix
U = zeros(numTx, numTx);
I = zeros(numTx, 1);
for l=1:numTx
    I(l) = (l-1) - (numTx-1)/2;
end
for i=1:numTx
    for j=1:numTx
        U(i,j) = exp(-1i*2*pi*I(j)/numTx*I(i))/sqrt(numTx);
    end
end

H_b = ctranspose(U)*H; % Beamspace channel matrix
s = (randn(N,1) + 1i*randn(N,1))/sqrt(2); % symbol vertor

% Compute beamspace channel power
for i=1:size(H_b, 2)
    H_b_power =conj(H_b).*H_b;
end

% contour(I, 1:N, transpose(H_b_power));
% plot3(I, 1:N, H_b_power);
% plot(conj(H_b).*H_b);
for i=1:N
    plot(conj(H_b(:,i)).*H_b(:,i))
end

% Precoding Matrix
for i=1:N
    F_MF = H; % Matched filter
    F_ZF = H/(ctranspose(H)*H); % Zero-forcing 
    eta = noise_power*numRx/Tx_power; % Regulization factor pf wiener filter
    F_WF = (H*ctranspose(H)+eta*eye(numRx))\H; % Wiener filter

    alpha_MF = sqrt(Tx_power/trace(F_MF*ctranspose(F_MF)));
    alpha_ZF = sqrt(Tx_power/trace(F_ZF*ctranspose(F_ZF)));
    alpha_WF = sqrt(Tx_power/trace(F_WF*ctranspose(F_WF)));

    G_MF = alpha_MF * F_MF;
    G_ZF = alpha_ZF * F_ZF;
    G_WF = alpha_WF * F_WF;
    
    % Beamspace channel matrix and precoding matrix
    
    G_b_MF = U\G_MF;
    G_b_ZF = U\G_ZF;
    G_b_WF = U\G_WF;
    
    r_MF = ctranspose(H)*G_MF*s;
    r_b_MF = ctranspose(H_b)*(G_b_MF)*s; %+ noise_power*(randn(N,1)+1i*randn(N,1))/sqrt(2);
end







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
