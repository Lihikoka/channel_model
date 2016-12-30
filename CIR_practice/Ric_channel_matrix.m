function H = Ric_channel_matrix(numTx, numRx, K_dB, L)
% Rician channel matrix
% Input : K_dB = K factor[dB]
% Output: H = 3-D Channel Matrix
K = 10^(K_dB/10);
H = sqrt(K/(K+1)) + sqrt(1/(K+1))*(randn(numTx, numRx,L)+1i*randn(numTx, numRx ,L))/sqrt(2);
