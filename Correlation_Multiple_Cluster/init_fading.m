function FadingMatrixTime = init_fading(FadingLength, ... 
					FadingOversamplingFactor, ... 
					NumberOfTxAntennas, ... 
					NumberOfRxAntennas, ... 
					NumberOfPaths, FadingType, sigma_rad, ... 
                    angle_rad) 
 
% FadingMatrixTime = init_fading(FadingLength, 
% FadingOversamplingFactor, NumberOfTxAntennas, NumberOfRxAntennas, 
% NumberOfPaths, FadingType, sigma_rad, angle_rad) 
% 
% Generates a matrix of fading coefficients of size 
% (NumberofTxAntennas*NumberOfRxAntennas*NumberOfPaths) x 
% FadingLength spanning FadingLength/FadingOversamplingFactor * 
% wavelength/2*speed seconds with a sampling time of 
% wavelength/(2*FadingOversamplingFactor*speed). The Doppler 
% spectrum of the fading coefficients can be either classic, 
% flat or "Laplacian". 
% 
% Inputs 
% 
% * Variable FadingLength, number of fading samples generated 
%   per tap coefficient 
% * Variable FadingOversamplingFactor, oversampling factor of 
%   the fading process in the time domain, so as to achieve 
%   a sampling time of wavelength/(2*FadingOversamplingFactor*speed). 
% * Variable NumberOfTxAntennas, number of antenna elements 
%   at Tx 
% * Variable NumberOfRxAntennas, number of antenna elements 
%   at Rx 
% * Variable NumberOfPaths, number of taps of the PDP 
% * Variable FadingType, type of the Doppler spectrum, which 
%   can be classic, flat or "Laplacian" 
% * Variable sigma_rad, standard deviation of the Laplacian 
%   PAS in the case of a "Laplacian" Doppler spectrum 
% * Variable angle_rad, AoA of the waves in the case of a 
%   Laplacian PAS 
% 
% Output 
% 
% * 2-D matrix FadingMatrix Time of size 
%   (NumberofTxAntennas*NumberOfRxAntennas*NumberOfPaths) x 
%   FadingLength containing FadingLength samples of 
%   NumberofTxAntennas * NumberOfRxAntennas * NumberOfPaths 
%   independent fading processes exhibiting classic, flat or 
%   "Laplacian" Doppler spectrum. 
% 
% 
% STANDARD DISCLAIMER 
% 
% CSys is furnishing this item "as is". CSys does not provide any 
% warranty of the item whatsoever, whether express, implied, or 
% statutory, including, but not limited to, any warranty of 
% merchantability or fitness for a particular purpose or any 
% warranty that the contents of the item will be error-free. 
% 
% In no respect shall CSys incur any liability for any damages, 
% including, but limited to, direct, indirect, special, or 
% consequential damages arising out of, resulting from, or any way 
% connected to the use of the item, whether or not based upon 
% warranty, contract, tort, or otherwise; whether or not injury was 
% sustained by persons or property or otherwise; and whether or not 
% loss was sustained from, or arose out of, the results of, the 
% item, or any services that may be provided by CSys. 
% 
% (c) Laurent Schumacher, AAU-TKN/IES/KOM/CPK/CSys - February 2002 
 
CutOff = floor(FadingLength/(2*FadingOversamplingFactor))-1; 
if strcmp(FadingType, 'laplacian') % Laplacian constrained PAS 
    FadingMatrixFreq = []; 
    for ii = 1:size(sigma_rad,2) 
        tmp = sqrt(1/(pi*FadingLength))... 
            ./sqrt((1+1e-9)... 
            .*ones(size(-CutOff:1:CutOff))-((-CutOff:1:CutOff)./CutOff).^2)... 
            .*(exp(-sqrt(2)/sigma_rad(ii).*abs(acos((-CutOff:1:CutOff)./CutOff)-angle_rad(ii)))... 
               + exp(-sqrt(2)/sigma_rad(ii).*abs(acos((-CutOff:1:CutOff)./CutOff)+angle_rad(ii)))); 
           FadingMatrixFreq = [FadingMatrixFreq; 
                               ones(NumberOfTxAntennas*NumberOfRxAntennas,1)... 
                               *[tmp(CutOff+2:2*CutOff-1),... 
                                 zeros(1,FadingLength-2*CutOff+3),... 
                                 tmp(3:CutOff+1)]]; 
    end; 
else 
    if strcmp(FadingType, 'classic') % Classical Doppler spectrum 
        tmp = sqrt(1/(pi*FadingLength))... 
            ./sqrt((1+1e-9)... 
            .*ones(size(1:1:CutOff))-((1:1:CutOff)./CutOff).^2); 
    elseif strcmp(FadingType, 'flat') % Flat Doppler spectrum 
        tmp = ones(size(1:1:CutOff)); 
    end; 
    FadingMatrixFreq = ones(NumberOfTxAntennas... 
        *NumberOfRxAntennas... 
		*NumberOfPaths,1)... 
        *[tmp(1:CutOff-1),... 
          zeros(1,FadingLength-2*CutOff+3),... 
          fliplr(tmp(2:CutOff-1))]; 
end; 
% Addition of a random phase 
FadingMatrixFreq = FadingMatrixFreq... 
                   .*exp((2*pi*i).*rand(NumberOfTxAntennas... 
					*NumberOfRxAntennas... 
					*NumberOfPaths, ... 
					FadingLength)); 
FadingMatrixTime = ifft(FadingMatrixFreq,FadingLength,2); 
% Normalisation to 1 
FadingMatrixTime = FadingMatrixTime./sqrt(mean(abs(FadingMatrixTime).^2,2)*... 
				      ones(1,size(FadingMatrixTime,2)));
