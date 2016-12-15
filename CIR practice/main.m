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

AoAMean = 0;
for i = 1:N
    AoAMean = AoAMean + mean(CIR(i).AOAs)/N;
end

AoAVariance = 0;
for i = 1:N
    AoAVariance = AoAVariance + 
end


