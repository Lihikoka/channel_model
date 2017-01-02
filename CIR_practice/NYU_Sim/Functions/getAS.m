function AS = getAS(aziArray,powerArray)

    deltaArray = (0:359)';
    ASMat_Temp = -ones(size(deltaArray));
    for deltaIdx = 1:numel(deltaArray)
        
        delta = deltaArray(deltaIdx);
        angle_temp = mod(aziArray+delta,360);
        
        meanAzi = sum(angle_temp.*powerArray) / sum(powerArray);
        mean_AziSq = sum(angle_temp.^2.*powerArray) / sum(powerArray);
        RMSAS = sqrt(mean_AziSq - meanAzi^2);
        
        if numel(powerArray) == 1
            RMSAS = 0;
        else
        end
        ASMat_Temp(deltaIdx) = RMSAS;
        
    end%end of deltaIdx
    
    AS = min(ASMat_Temp);
    
end