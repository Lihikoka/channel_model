clc
clear
load('Results_2873_GHz_Outdoor_20_CIRs.mat')
N = 20;
sum = 0;
tap = 2;
HPowerSum = 0;
for i = 1:81
    HPowerSum = HPowerSum + CIR_Struct.CIR_1.HPowers{tap, 1}(i);
end
CIR_Struct.CIR_1.pathPowers(tap)
HPowerSum
HPowerSum - CIR_Struct.CIR_1.pathPowers(tap)