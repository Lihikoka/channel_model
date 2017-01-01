clc
clear

n = 4;
K = 100;
a = zeros(4, 1);
I = zeros(n, 1);
for l=1:n
    I(l) = (l-1) - (n-1)/2;
end