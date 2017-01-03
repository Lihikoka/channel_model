function CIR_Struct = getLocalCIR(CIR_Struct,TxArrayType,RxArrayType,Nt,Nr,Wt,Wr,dTxAnt,dRxAnt) 
% This function generates the local area CIRs. 
%
% Inputs:
%   - TxArrayType: a string, specifying the type of TX antenna array 
%       (e.g., 'ULA', 'URA')
%   - RxArrayType: a string, specifying the type of RX antenna array
%       (e.g., 'ULA', 'URA')
%   - Nt: the number of TX antennas
%   - Nr: the number of RX antennas 
%   - Wt: the number of TX antennas in the azimuth dimention
%   - Wr: the number of RX antennas in the azimuth dimention
%   - dTxAnt: the spacing between adjacent TX antennas in units of wavelengths
%   - dRxAnt: the spacing between adjacent TX antennas in units of wavelengths 
% Output:
%   - CIR: structure containing the MIMO channel coefficients

% Example input parameters (for test purpose only): 
% TxArrayType = 'ULA'; RxArrayType = 'URA'; Nt = 4; Nr = 10; ...
%   dTxAnt = 1/2; dRxAnt = 1/2; Wt = 1; Wr = 2;
%
% Copyright ?2016 NYU

% imaginary unit
j = sqrt(-1);

% speed of light (m/s)
c = 3e8; 

% Rician K factor
K = 10^(10/10); 

% Number of CIRs to process
N = size(fieldnames(CIR_Struct),1);
for CIRIdx = 1:N

    % Extract n-th CIR
    CIR = CIR_Struct.(['CIR_',num2str(CIRIdx)]);
    
    % carrier frequency (GHz)
    f = str2double(CIR.frequency(1:2))*1e9; 

    % carrier wavelength
    wl = c/f;    
    
    % determine if the environment is LOS
    if strcmp(CIR.environment,'LOS') == true 
        % spatial correlation coefficients for LOS
        A = 0.99; B = 1.95; C = 0; 
    elseif strcmp(CIR.environment,'NLOS') == true
        % spatial correlation coefficients for NLOS
        A = 0.9; B = 1; C = -0.1; 
    else
    end

    % number of paths
    nTap = length(CIR.pathDelays); 

    % absolute propagation time delay of each path
    timeDelay = CIR.pathDelays.*1e-9; 

    % received power in mW of each path
    Pr_lin = CIR.pathPowers; 

    % Calculate local area CIRs and the corresponding parameters in a path-by-path manner
    Hw_RC = cell(nTap,1); % Yen-hua modified
    H_RC = cell(nTap,1); % Yen-hua modified
    for a = 1:nTap
        % Determine the Tx and Rx antenna array types
        if strcmp(TxArrayType,'ULA') == true && strcmp(RxArrayType,'ULA') == true 

            % Calculate Rx spatial correlation matrix Rr
            for i = 1:Nr  
                for m = 1:Nr
                    % Spatial correlation between Rx antenna elements i and m in a ULA
                    Rr(i,m) = (A.*exp(-B.*dRxAnt.*wl.*abs(i-m))-C).*exp(-j.*unifrnd(-pi,pi).*(i-m)); 
                end
            end 

            % Calculate Tx spatial correlation matrix Rt
            for i = 1:Nt  
                for m = 1:Nt
                    % Spatial correlation between Tx antenna elements i and m in a ULA
                   Rt(i,m) = (A.*exp(-B.*dTxAnt.*wl.*abs(i-m))-C).*exp(-j.*unifrnd(-pi,pi).*(i-m)); 
                 end
            end 

            % Determine the Tx and Rx antenna array types
            elseif strcmp(TxArrayType,'ULA') == true && strcmp(RxArrayType,'URA') == true 
                % Calculate Rx spatial correlation matrix Rr
                for i = 1:Nr  
                    for m = 1:Nr
                       % Spatial correlation between Rx antenna elements i and m in a URA
                       Rr(i,m) = (A.*exp(-B.*dRxAnt.*wl.*sqrt((mod(i,Wr)-mod(m,Wr))^2+(fix(i/Wr)-fix(m/Wr))^2))-C).*exp(-j.*unifrnd(-pi,pi).*(i-m)); 
                    end
                end

            % Calculate Tx spatial correlation matrix Rt
            for i = 1:Nt  
                for m = 1:Nt
                   % Spatial correlation between Tx antenna elements i and m in a ULA
                   Rt(i,m) = (A.*exp(-B.*dTxAnt.*wl.*abs(i-m))-C).*exp(-j.*unifrnd(-pi,pi).*(i-m)); 
                end
            end

            % Determine the Tx and Rx antenna array types
            elseif strcmp(TxArrayType,'URA') == true && strcmp(RxArrayType,'ULA') == true
            % Calculate Rx spatial correlation matrix Rr
                for i = 1:Nr  
                    for m = 1:Nr
                       % Spatial correlation between Rx antenna elements i and m in a ULA
                       Rr(i,m) = (A.*exp(-B.*dRxAnt.*wl.*abs(i-m))-C).*exp(-j.*unifrnd(-pi,pi).*(i-m)); 
                    end
                end

            % Calculate Tx spatial correlation matrix Rt
            for i = 1:Nt  
                for m = 1:Nt
                   % Spatial correlation between Tx antenna elements i and m in a URA
                   Rt(i,m) = (A.*exp(-B.*dTxAnt.*wl.*sqrt((mod(i,Wt)-mod(m,Wt))^2+(fix(i/Wt)-fix(m/Wt))^2))-C).*exp(-j.*unifrnd(-pi,pi).*(i-m)); 
                 end
            end

            % Determine the Tx and Rx antenna array types
            elseif strcmp(TxArrayType,'URA') == true && strcmp(RxArrayType,'URA') == true
            % Calculate Rx spatial correlation matrix Rr
                for i = 1:Nr  
                for m = 1:Nr
                    % Spatial correlation between Rx antenna elements i and m in a URA
                   Rr(i,m) = (A.*exp(-B.*dRxAnt.*wl.*sqrt((mod(i,Wr)-mod(m,Wr))^2+(fix(i/Wr)-fix(m/Wr))^2))-C).*exp(-j.*unifrnd(-pi,pi).*(i-m)); 
                 end
                end

            % Calculate Tx spatial correlation matrix Rt
            for i = 1:Nt  
                for m = 1:Nt
                   % Spatial correlation between Tx antenna elements i and m in a URA
                   Rt(i,m) = (A.*exp(-B.*dTxAnt.*wl.*sqrt((mod(i,Wt)-mod(m,Wt))^2+(fix(i/Wt)-fix(m/Wt))^2))-C).*exp(-j.*unifrnd(-pi,pi).*(i-m)); 
                 end
            end

        end%% end of if statement

        % eigenvalue decomposition of Rr and Rt
        [Ur,Dr] = eig(Rr); [Ut,Dt] = eig(Rt); 

        % random Rayleigh distribution (to be incoperated in the Rician distribution)
        randRayleigh = randn(Nr,Nt)+j*randn(Nr,Nt); 

        % H matrix with independent Rician distribution
        Hw_RC{a,1} = sqrt(K/(K+1))*exp(j*pi/4)+sqrt(1/2/(K+1))*randRayleigh; 
        
        % H matrix for the a-th path, where the elements obey the small-scale Rician distribution
        % specified by Hw_RC and spatial correlation specified by Rt and Rr
        H_RC{a,1} = sqrt(Pr_lin(a))*Ur*Dr^(1/2)*Hw_RC{a,1}*Dt^(1/2)*Ut'; 
        
        % time delay for the a-th path
        CIR.HDelays{a,1} = CIR.pathDelays(a); 

        % received power between each Tx antenna and each Rx antenna for the a-th path
        CIR.HPowers{a,1} = abs(H_RC{a,1}).^2; 
        
        % phase between each Tx antenna and each Rx antenna for the a-th path
        CIR.HPhases{a,1} = CIR.pathPhases(a)+angle(H_RC{a,1}); 
        CIR.HAODs{a,1} = CIR.AODs(a); % AOD for the ath path
        CIR.HZODs{a,1} = CIR.ZODs(a); % ZOD for the ath path
        CIR.HAOAs{a,1} = CIR.AOAs(a); % AOA for the ath path
        CIR.HZOAs{a,1} = CIR.ZOAs(a); % ZOA for the ath path

    end% end of aTap for loop

    % number of antenna elements in the Tx array
    CIR.NumOfTxElements = Nt; 

    % number of antenna elements in the Rx array
    CIR.NumOfRxElements = Nr; 

    % H matrix
    CIR.H = H_RC; 

    % Store the updated CIR
    CIR_Struct.(['CIR_',num2str(CIRIdx)]) = CIR;
    
end%end of CIRIdx
whos
end% end of function

