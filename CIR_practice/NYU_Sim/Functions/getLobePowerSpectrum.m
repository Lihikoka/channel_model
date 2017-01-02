function lobePowerSpectrum_struct = ...
    getLobePowerSpectrum(numberOfLobes,cluster_subpath_lobe_mapping,powerSpectrum,angleType)
% Generate the lobe power angular spectra
%
% Inputs:
%   - numberOfLobes: the number of spatial lobes
%   - cluster_subpath_lobe_mapping: a structure containing the mapping
%   between cluster number, subpath number, and spatial lobe number for
%   each multipath 
%   - powerSpectrum: an array containing all multipath parameters
%   - angleType: 'AOD' or 'AOA'
% Output:
%   - lobePowerSpectrum_struct: a structure containing lobe angular spectra
%
% Copyright © 2016 NYU

if strcmp(angleType,'AOD') == true
    subpathSpectrum = powerSpectrum(:,1:5);
elseif strcmp(angleType,'AOA') == true
    subpathSpectrum = powerSpectrum(:,[1:3 6:7]);
else
end

%%% initialize lobe power spectrum
lobePowerSpectrum_struct = struct;

for lobeIndex = 1:numberOfLobes
    
   indSameLobe = find(cluster_subpath_lobe_mapping(:,3) == lobeIndex);
   
   %%% subpaths that belong to same lobe
   subpathSpectrum_SameLobe = subpathSpectrum(indSameLobe,:);
    
   %%% store
   lobePowerSpectrum_struct.(['Lobe',num2str(lobeIndex)]) = subpathSpectrum_SameLobe;
    
end



end