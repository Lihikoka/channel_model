clc
clear

frequency = '2873_GHz';
scenario = 'outdoor';
N=100;
K_dB = 5; % 5~15

string = ['load Results_',frequency,'_',scenario,'_',num2str(N),'_CIRs'];
eval(string);

allAOA = [];
allPower = [];
indexLOS = [];




%%%%%%%%% fetch CIR data
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
end

% Generate Uncorrelated Rician Fading channel matrix
for i=1:1
    L = length(CIR(indexLOS(i)).pathDelays);
    A = Ric_channel_matrix(numTx, numRx, K_dB, L);
end

% Generate Tx and Rx correlation matrix
Rtx = zeros(numTx, numTx);
Rrx = zeros(numRx, numRx);
for i=1:numRx
    for j=i:numRx
        %Rrx(i, j) = 
    end
end


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
