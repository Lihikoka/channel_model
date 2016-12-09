clc
clear

frequency = '28_GHz';
scenario = 'outdoor';
N = 100;
cluster_number = 1;

string = ['load Results_',frequency,'_',scenario,'_',num2str(N),'_CIRs'];
eval(string);

CIR = CIR_Struct.CIR_1;
fields = fieldnames(CIR);
hh = CIR.('H')(1);
hh{1, 1};
P_l = CIR.('pathPowers');

A_l = randn(4);
a_l = reshape(A_l, [1, size(A_l, 1)^2]);
