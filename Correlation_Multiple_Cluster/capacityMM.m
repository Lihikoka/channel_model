function sum_inst_cap = capacityMM(nTx, nRx, H, e, V, Plinear, sigma_N_2); 
% sum_inst_cap = capacity(nTx, nRx, H, e, s, V, Plinear, sigma_N_2) 
% 
% Computes the capacity of the channel from the EVD of the H matrix 
% and the SVD of the correlation matrix R using Uniform power allocation 
 
Le      = min(nTx, nRx); % ominaisarvojen max-m??r? (Rank(R)) 
Signal  = zeros(nTx, nTx, size(e,2)); 
Noise   = sigma_N_2*eye(nRx);    
 
% 
% Computation of instantaneous capacities 
% 
 
% Brute force method 
% C = det ( I + H Q H') where Q = V Pw V' and H'H = V D V' 
% Note: Telatar writes H'H = U' D U and Q = U' Pw U 
 
for ii = 1:size(H,3) 
  icu(ii) = log2(det(eye(nRx) + inv(Noise) * H(:,:,ii) * ((Plinear/nTx) ... 
			* eye(nTx)) * H(:,:,ii)')); 
end; 
 
% 
% Computation of cumulative capacities 
% 
 
sum_inst_cap(1,:) = mean(real(icu)); % Uniform capacities (Brute-Force method) 