clc
clear

n = 81;
K = 20;
a = zeros(n, n);
I = zeros(n, 1);
for l=1:n
    I(l) = (l-1) - (n-1)/2;
end


for i=1:n
    for j=1:n
        a(i,j) = exp(-1i*2*pi*I(j)/n*I(i)) /sqrt(n);
    end
end