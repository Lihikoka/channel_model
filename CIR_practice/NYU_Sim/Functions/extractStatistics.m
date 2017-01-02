function extractStatistics(CIR_Struct)
% This function extract second-order time and angular statistics, and saves
% them in a .txt file in the \Results folder.
%
% Inputs:
%   - CIR_Struct: a structure containing N independently generated CIRs
% Output:
%   - A .txt file containing the RMS delay spreads, ZSD/ASD/ZSA/ASA values.
%   .txt file save in the \Results folder.
%
% Copyright © 2016 NYU

% Output to command window
disp('Extracting 2nd order statistics...')

% Parameters needed for filename of output .txt file
frequency = CIR_Struct.CIR_1.frequency;
scenario = CIR_Struct.CIR_1.scenario;
N = size(fieldnames(CIR_Struct),1);

% Create the output file name
saveName = ['SecondOrderStats_',frequency,'_',scenario,'_',num2str(N),'_CIRs.txt'];

% Open the file identifier
fid = fopen(saveName,'w');

% Print first line
fprintf(fid,'%s\t%s\t%s\t%s\t%s\t%s\t%s\n','CIR#','Environment','RMS DS (ns)','ZSD (deg)','ASD (deg)','ZSA (deg)','ASA (deg)');

for CIRIdx = 1:N
    
    % Extract multipath parameters
    pathDelays = CIR_Struct.(['CIR_',num2str(CIRIdx)]).pathDelays;
    pathPowers = CIR_Struct.(['CIR_',num2str(CIRIdx)]).pathPowers;
    AODs = CIR_Struct.(['CIR_',num2str(CIRIdx)]).AODs;
    ZODs = CIR_Struct.(['CIR_',num2str(CIRIdx)]).ZODs;
    AOAs = CIR_Struct.(['CIR_',num2str(CIRIdx)]).AOAs;
    ZOAs = CIR_Struct.(['CIR_',num2str(CIRIdx)]).ZOAs;
    
    % Compute RMS DS
    meanTau = sum(pathDelays.*pathPowers)/sum(pathPowers);
    mean_TauSq = sum(pathDelays.^2.*pathPowers)/sum(pathPowers);
    RMSDS = sqrt(mean_TauSq - meanTau^2);
    if numel(pathDelays) == 1
        RMSDS = 0;
    else
    end

    % Compute AOD/AOA Angular Spreads
    ASD = getAS(AODs,pathPowers);
    ASA = getAS(AOAs,pathPowers);
    
    % Compute ZOD
    ZSD = getAS(ZODs,pathPowers);
    ZSA = getAS(ZOAs,pathPowers);
    
    % Extract environment
    environment = CIR_Struct.(['CIR_',num2str(CIRIdx)]).environment;
    
    % print to .txt
    fprintf(fid,'%s\t%s\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n',['CIR',num2str(CIRIdx)],environment,RMSDS,ZSD,ASD,ZSA,ASA);
    
end

% Close file identifier
fclose('all');

% Output to command window
disp('Extracting 2nd order statistics...Done.')
disp('Statistics saved.')



end