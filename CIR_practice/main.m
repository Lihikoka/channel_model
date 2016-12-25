clc
clear

frequency = '28_GHz';
scenario = 'outdoor';
N=100;
string = ['load Results_',frequency,'_',scenario,'_',num2str(N),'_CIRs'];
eval(string);

allAOA = [];
for i = 1:N
    string = ['CIR(', num2str(i), ') = ', 'CIR_Struct.CIR_', num2str(i), ';'];
    eval(string);
    allAOA = [allAOA; CIR(i).AOAs];
end

AOAmean = mean(allAOA);
AOAStdev = std(allAOA);
% delta_phi_deg = AOAStdev * sqrt(3);

%%%%%%%%%% plot AoA distribution
AOARange = 0:0.01:360;
AOAGraph = zeros(1, 36001);

for i = 1:length(allAOA)
    pos          = find(AOARange>=allAOA(i)); 
    AOAIndex = ((AOARange(pos(1))-allAOA(i))*(pos(1)-1) + ... % interpolation
                (allAOA(i)-AOARange(pos(1)-1))*(pos(1))) / ...
                (AOARange(pos(1))-AOARange(pos(1)-1));
    AOAIndex = round(AOAIndex);
    AOAGraph(AOAIndex) = AOAGraph(AOAIndex) + 1;
end

figure(1)
plot(AOARange, AOAGraph)
xlabel('AOA(deg)');
ylabel('number');

