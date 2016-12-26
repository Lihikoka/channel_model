clc
clear

frequency = '73_GHz';
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
AOApdf = zeros(1, 36001);

for i = 1:length(allAOA)
    pos          = find(AOARange>=allAOA(i)); 
%     AOAIndex = round(AOARange(pos(1)));
    AOAIndex = ((AOARange(pos(1))-allAOA(i))*(pos(1)-1) + ... % interpolation
                (allAOA(i)-AOARange(pos(1)-1))*(pos(1))) / ...
                (AOARange(pos(1))-AOARange(pos(1)-1));
    AOAIndex = round(AOAIndex);
    AOAGraph(AOAIndex) = AOAGraph(AOAIndex) + 1;
end

sum = 0;
for i = 1:length(AOARange)
    sum = sum + AOAGraph(i);
    if ~mod(AOARange(i), 1)
        AOApdf(i) = sum/length(allAOA);
        sum=0;
    end
end


figure(1)
plot(AOARange, AOApdf)
xlabel('AOA(deg)');
ylabel('number');

