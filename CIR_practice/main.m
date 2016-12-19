clc
clear

frequency = '28_GHz';
scenario = 'outdoor';
N=100;
string = ['load Results_',frequency,'_',scenario,'_',num2str(N),'_CIRs'];
eval(string);

for i = 1:N
    string = ['CIR(', num2str(i), ') = ', 'CIR_Struct.CIR_', num2str(i), ';'];
    eval(string);
end

AoASum = 0;
pathAmount = 0;
for i = 1:N
    AoASum = AoASum + sum(CIR(i).AOAs);
    pathAmount = pathAmount + length(CIR(i).AOAs);
end
AoAMean = AoASum / pathAmount;

AoAVariance = 0;
for i = 1:N
    AoAVariance = AoAVariance + sum((CIR(i).AOAs - AoAMean).^2);
end
AoAStd = sqrt(AoAVariance/pathAmount);
delta_phi_deg = AoAStd * sqrt(3);

%%%%%%%%%%% plot AoA distribution
AoArange = -180:0.01:180;
AoAGraph = zeros(1, 36001);
for i = 1:N
    for ii = 1:36001
        AoA_index = ((tmp(pos(1),1)-AS_deg(k))*tmp(pos(1)-1,2) + ...
                    (AS_deg(k) - tmp(pos(1)-1,1))*tmp(pos(1),2))/ ...
                    (tmp(pos(1),1)-tmp(pos(1)-1,1));
        AoA_index = (CIR(i).AOAs-AoArange(ii))
    end
end



